# cpp_Standard

The full C++17 Fourier_Quad pipeline build (includes PCA `PSFRecons`).

The complete C++ pipeline guide - source structure, pipeline stages,
configuration, building, run modes, initializer output layout, Docker, and HPC -
now lives in [`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference
is [`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

## Runtime configuration

`RuntimeConfig` now owns every run-selectable Process, ExtCat, Init, and
Standard Lensing value. Start from [`pipeline.example.ini`](pipeline.example.ini)
and launch with:

```bash
./Fourier_Quad_Pipe --config pipeline.ini
```

The precedence is compiled header defaults, then the INI file, then CLI
overrides. Both `--name value` and `--name=value` remain supported. Rank 0 alone
reads the file and broadcasts its exact text; every MPI rank parses and validates
the same text before the write-once configuration store is initialized.

The INI parser accepts `[process]`, `[extcat]`, `[init]`, and `[lensing]`, `#` or
`;` comments, quoted strings, comma-separated lists, and common boolean forms.
Unknown sections/keys and malformed values fail startup with file, section, key,
value, and reason. `source_cat` is a compatibility alias for the authoritative
`extcat.output_directory` field.

`ext_cat=0` selects the internal Stage-3/9 catalog path and its 30-column rearr
schema without resolving or configuring an external catalog. FD still requires
`ext_cat=1` and fails validation otherwise. `nmax_chip`, `chipnx`, and `chipny`
are runtime values with positive-value validation but no DECam-specific hard
ceiling. The removed `npx`, `npy`, and `NMAX_EXPO` settings are not accepted.

## Stamp-cube format

Stage 3--7 stamp collections are contiguous three-dimensional FITS images with
axes `(x, y, stamp)`. Their in-memory layout is `[stamp][row][col]`, with flat
index `((stamp * ny) + row) * nx + col`. Writers add
`FQFMT='STAMP_CUBE'` and `FQORDER='X,Y,STAMP'`; readers recover all dimensions
from the FITS header, and each pipeline consumer checks them against its catalog
count and configured stamp size.

This is an intentional format break. Products written as legacy two-dimensional
stamp mosaics are rejected, and Standard PCA products generated before this
row-major cleanup must be regenerated from Stage 3 onward.

## Norm FITS background and sigma metadata

Stage 1 writes the final background and sigma-plane coefficients into the
existing `*_norm.fits` primary header together with the normalized pixels. For
`ccd_split=1`, the logical long-string keywords are `BGCO` and `SIGCO`; for
`ccd_split=2`, they are `BG1CO`, `BG2CO`, `SIG1CO`, and `SIG2CO`. Coefficients
are serialized with scientific precision 17, and Stage 3 reads the image and
all coefficient keywords through one FITS handle.

Stage 3 reconstructs the sigma plane from the header and subtracts the Stage-1
background model once per amplifier after the existing optional flat correction.
The old sigma metadata pixels are no longer read or written; missing or
malformed header metadata is a chip failure with no legacy-pixel fallback.

## Build and focused verification

Use the portable Makefile inputs described in the main guide. For example, when
MPI/scientific libraries share one prefix but Eigen uses another include path:

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
  EIGEN_INCLUDE="${EIGEN_INCLUDE}" test-psf-recons-orientation
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
and metadata, and dense selected-plane ordering. The PSF test exercises the same
row-major copy, centering, and reconstruction helpers used by `PSFRecons`, using
asymmetric pixels and a numerical comparison with the legacy feature permutation.
The runtime-config test covers compiled defaults, every INI value family,
transactional diagnostics, file/CLI precedence, list behavior, stage and FD
dependencies, non-default geometry, internal rearr schema, and the write-once
store. The catalog-row-count test covers equal 100-row catalogs, both mismatch
directions, and paired header-only catalogs. The PSF-state test covers actual
candidate counts of 0, 10, 300, and 2301, including full live-stride matrix
access beyond the 2000-row reservation hint.

The current local WSL2 verification stack is GCC/G++ 15.2.0, Open MPI 5.0.10,
CFITSIO 4.6.3, FFTW 3.3.10, Eigen 3.4.0, and OpenBLAS/LAPACK 0.3.33. The portable
cluster target and module-compatible versions are recorded in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md#compiler-and-libraries).

## Dynamic catalog and PSF capacities

`LensingConfig::ngal_max` (4000) and `LensingConfig::nstar_max` (2000) are
initial metadata-vector reservation capacities, not hard catalog limits. Stage
3 appends every accepted source and star candidate, while the large stamp
buffers continue to grow on demand rather than reserving the full hint size.

Stage 5 creates empty state for the exposure's live chips, loads each chip's
complete candidate catalog, and then allocates a full
`actual_nstar x actual_nstar` pairwise chi matrix. The pairwise calculation,
star selection, and local/hybrid PSF fitting algorithms are unchanged. In the
external-catalog Stage-9 branch, the physical data-row counts after the two
headers must match before combination; a mismatch aborts the MPI world, while
the existing row parsing, cuts, calibration, and pairing order remain unchanged.

## Runtime catalog layout

External `CatalogLayout` is resolved once in `main` from `RuntimeConfig` and
passed only to phases that consume external rows. Rearrangement receives the
smaller external-or-internal `RearrCatalogSchema`. In
pass-through mode, `EXTCAT_TOTAL_COLUMNS` supplies the external prefix width.
With explicit projection, the external width is the projection length and all
downstream CCD, source-suffix, and complete-row offsets follow that runtime
width. The fixed process-main suffix remains 29 fields.

Explicit projections must preserve RA, Dec, and ZP. Each g/r/i/z/y magnitude
is optional: raw configured column `0`, or omission of a positive configured
identity from a projection, marks that band absent. Other physical extcat
columns may remain in a row and contribute to its width, but have no
`CatalogLayout` member and are not consumed by downstream algorithms. FD
selects one available magnitude in i -> z -> r -> g -> y order before catalog
I/O and fails if none exists; the selected band drives both the magnitude-range
cut and the size-magnitude star bar. The focused `test-catalog-layout` target
covers the legacy 48-column schema, a minimal 3 + 1 + 29 = 33-column schema,
optional/reordered projections, all FD fallbacks, selected-band ingestion,
extra unmodeled fields, invalid projections, and exact FD row-width validation.

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
