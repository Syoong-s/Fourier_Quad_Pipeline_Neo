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
