# cpp_Lite

The frozen-branch simplified C++17 Fourier_Quad pipeline build (`PSFRecons`
removed). See `REFACTOR_NOTES.md` for the Lite change log.

The complete C++ pipeline guide lives in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md). The full parameter reference is
[`../CPP_PIPELINE_PARAMETERS.md`](../CPP_PIPELINE_PARAMETERS.md).

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
./Fourier_Quad_Pipe --help
```

The stamp-cube test verifies a non-square `5 x 3 x 4` round trip, raw FITS axes
and metadata, and dense selected-plane ordering. The current local WSL2 stack is
GCC/G++ 15.2.0, Open MPI 5.0.10, CFITSIO 4.6.3, FFTW 3.3.10, Eigen 3.4.0, and
OpenBLAS/LAPACK 0.3.33; portable cluster versions are recorded in
[`../CPP_GUIDE.md`](../CPP_GUIDE.md#compiler-and-libraries).
