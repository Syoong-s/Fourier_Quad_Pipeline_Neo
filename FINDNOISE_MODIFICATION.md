# Stage-3 `findNoise()` noise-candidate update

This change removes maximum-pixel ranking from Standard and Lite noise-stamp selection. The
existing 16-candidate perimeter geometry and chip-edge cut remain in place. Candidates now pass
fixed amplifier, mask, local-sigma, MAD, positive-tail, and final masked-fraction gates, are
shuffled with the existing rank-local `NumericalRecipes::ran1()` stream, and are accepted in that
random order. `markNoise()` and `decorateStamp()` use the accepted candidate's local sigma.

The same-amplifier gate uses the effective runtime `lensing.ccd_split` value. This keeps candidate
selection aligned with the amplifier geometry used to generate the chip's local-sigma map when an
INI file overrides the compiled fallback.

## Build environment

- Language standard: C++17.
- Local compiler: GCC/mpicxx 15.2.0 with Open MPI 5.0.10.
- Local libraries: CFITSIO 4.6.4, FFTW 3.3.11, OpenBLAS BLAS/LAPACK 3.11.0, and Eigen 3.4.0.
- Cluster target: a Linux MPI toolchain with a C++17 compiler and compatible CFITSIO, FFTW,
  BLAS/LAPACK, and Eigen development modules. Site-specific module names should be supplied by
  the cluster rather than hard-coded in this repository.

## Portable build and run

From either `cpp_Standard` or `cpp_Lite`:

```bash
make clean
make -j
./Fourier_Quad_Pipe --help
mpirun -np <ranks> ./Fourier_Quad_Pipe [options] [LEGACY_EXPO_LIST]
```

If Eigen is outside the compiler's default include path, pass its portable installation prefix:

```bash
make -j EIGEN_INCLUDE="$EIGEN_INCLUDE"
```

## Local verification override

The 2026-08-20 WSL2 verification used the repository Makefiles with the existing Pixi compiler
and library stack, plus `EIGEN_INCLUDE=/usr/include/eigen3`. Both variants completed a clean full
build and displayed their command-line help successfully.
