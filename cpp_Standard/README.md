# cpp_Standard

The full C++17 Fourier_Quad pipeline build (includes PCA `PSFRecons`).

The complete C++ pipeline guide - source structure, pipeline stages,
configuration, building, run modes, initializer output layout, Docker, and HPC -
now lives in [`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference
is [`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

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
./Fourier_Quad_Pipe --help
```

The stamp-cube test verifies a non-square `5 x 3 x 4` round trip, raw FITS axes
and metadata, and dense selected-plane ordering. The PSF test exercises the same
row-major copy, centering, and reconstruction helpers used by `PSFRecons`, using
asymmetric pixels and a numerical comparison with the legacy feature permutation.

The current local WSL2 verification stack is GCC/G++ 15.2.0, Open MPI 5.0.10,
CFITSIO 4.6.3, FFTW 3.3.10, Eigen 3.4.0, and OpenBLAS/LAPACK 0.3.33. The portable
cluster target and module-compatible versions are recorded in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md#compiler-and-libraries).

## Runtime catalog layout

`CatalogLayout` is resolved once in `main` from `RuntimeOptions` and passed to
`process_extcat`, `process_main`, `process_rearr`, and `process_fd`. In
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
