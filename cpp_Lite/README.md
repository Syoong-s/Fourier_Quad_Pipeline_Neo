# cpp_Lite

The frozen-branch simplified C++17 Fourier_Quad pipeline build (`PSFRecons`
removed). See `REFACTOR_NOTES.md` for the Lite change log.

The complete C++ pipeline guide lives in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference is
[`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

## Runtime configuration

`RuntimeConfig` owns every run-selectable Process, ExtCat, Init, and surviving
Lite Lensing value. Start from [`pipeline.example.ini`](pipeline.example.ini):

```bash
./Fourier_Quad_Pipe --config pipeline.ini
```

Precedence is compiled defaults, then INI values, then CLI overrides. Rank 0
alone reads the file and broadcasts its exact text; every MPI rank parses and
validates that text before the write-once store is initialized. Unknown keys,
sections, or malformed values fail startup transactionally with structured
diagnostics. `source_cat` aliases the authoritative
`extcat.output_directory`.

Lite exposes only `process_stage`, `astrometry_cat`, `ccd_split`, `gal_smooth`,
`star_smooth`, `pixel_size`, `nmax_chip`, `chipnx`, and `chipny` in
`[lensing]`. Its frozen Gaia/DQ/external-catalog/frame-star/deblending/local-PSF
behaviors remain implemented directly; Standard-only selector and path keys are
rejected instead of silently pretending to work. The removed `npx`, `npy`, and
`NMAX_EXPO` settings are not accepted.

## Stamp-cube format

Stage 3--7 stamp collections are contiguous three-dimensional FITS images with
axes `(x, y, stamp)`. Their in-memory layout is `[stamp][row][col]`, with flat
index `((stamp * ny) + row) * nx + col`. Writers add
`FQFMT='STAMP_CUBE'` and `FQORDER='X,Y,STAMP'`; readers derive all dimensions
from the FITS header, and consumers verify them against the configured stamp
size and the matching catalog count.

This intentionally breaks compatibility with legacy two-dimensional stamp
mosaics. Regenerate intermediate stamp products from Stage 3 onward. Lite keeps
its frozen `PSF_Ms=0` behavior and does not add Standard's PCA reconstruction.

## Build and focused verification

Use the portable Makefile inputs described in the main guide. For example:

```bash
make clean
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" -j4
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-stamp-cube-io
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-point-source-statistics
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-universalblock
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-catalog-layout
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-catalog-row-count
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-psf-model-state
make CXX="${MPI_PREFIX}/bin/mpicxx" \
  STACK_PREFIX="${STACK_PREFIX}" \
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-runtime-config
./Fourier_Quad_Pipe --help
```

The stamp-cube test verifies a non-square `5 x 3 x 4` round trip, raw FITS axes
and metadata, and dense selected-plane ordering. The catalog-row-count test
covers equal 100-row catalogs, both mismatch directions, and paired header-only
catalogs. The PSF-state test covers actual candidate counts of 0, 10, 300, and
2301, including full live-stride matrix access beyond the 2000-row reservation
hint.

The runtime-config test covers Lite defaults, all four INI sections, strict
rejection of deleted-branch keys, default/file/CLI precedence, geometry and
stage validation, transactional failures, and the immutable store.

The current local WSL2 stack is GCC/G++ 15.2.0, Open MPI 5.0.10, CFITSIO 4.6.3,
FFTW 3.3.10, Eigen 3.4.0, and OpenBLAS/LAPACK 0.3.33; portable cluster versions
are recorded in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md#compiler-and-libraries).

## Dynamic catalog and PSF capacities

`LensingConfig::ngal_max` (4000) and `LensingConfig::nstar_max` (2000) are
initial metadata-vector reservation capacities, not hard catalog limits. Stage
3 appends every accepted source and star candidate, while the large stamp
buffers continue to grow on demand rather than reserving the full hint size.

Stage 5 creates empty state for the exposure's live chips, loads each chip's
complete candidate catalog, and then allocates a full
`actual_nstar x actual_nstar` pairwise chi matrix. The pairwise calculation,
star selection, and local polynomial PSF fitting algorithms are unchanged.
Stage 9 requires the physical data-row counts after the shear and external
catalog headers to match before combination; a mismatch aborts the MPI world,
while the existing row parsing, cuts, calibration, and pairing order remain
unchanged.

## Runtime catalog layout

Lite uses the same startup-resolved external `CatalogLayout` contract as
Standard. `main` passes one immutable layout to `process_main`,
`process_rearr`, and `process_fd`; `process_extcat` consumes the same runtime
extcat section directly. Pass-through mode takes the external
prefix width from `EXTCAT_TOTAL_COLUMNS`; explicit projection uses the
projection length. CCD, source-suffix, and complete-row offsets are derived
from that effective width, while the process-main suffix remains 29 fields.

RA, Dec, and ZP are mandatory. Each g/r/i/z/y magnitude is optional: raw
configured column `0`, or omission of a positive configured identity from an
explicit projection, marks the band absent. Non-required physical extcat
columns contribute only to row width and have no layout member. At FD startup,
Lite selects i -> z -> r -> g -> y; the selected band drives the magnitude cut
and size-magnitude star bar, and FD fails before catalog I/O when no band is
available. `test-catalog-layout` covers layout mapping, fallback selection,
selected-band ingestion, and exact runtime row widths.

## Invalid norm chip gate

`namespace Universalblock` classifies the first pixel of each sharded
`*_norm.fits` as `Valid`, `Invalid`, `Missing`, or `ReadError`. Stages 3, 4, 5,
6, 7, and 9 check this status before opening their chip-level downstream inputs.
Only `Invalid` is a silent skip; missing or malformed norm FITS files remain
reported input errors, and valid chips retain their existing downstream-file
error paths.

Stage 5 keeps its established zero-star PSF placeholder behavior after an
invalid chip is skipped during candidate loading. Stage 9 leaves external-catalog
header discovery file-driven and applies the norm gate only in the data loop,
before opening `_shear.dat` or `_orig.cat`. The focused `test-universalblock`
target covers path derivation, valid and invalid sentinels, invalid-path silence,
missing norm files, and malformed FITS input.
