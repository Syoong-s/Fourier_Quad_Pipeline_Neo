# Star-candidate corrected-power quality gate

## Behavior

Stage 5 now validates every corrected star-candidate Fourier-power stamp before
shape diagnostics, pairwise chi construction, star clustering, selected-star
output, or PSF fitting.

A candidate remains active only when:

1. every power pixel is finite;
2. the median of the eight pixels surrounding the Fourier DC pixel is
   non-negative;
3. the signed total power is finite and strictly positive; and
4. the derived size, ellipticity, and FWHM diagnostics are finite.

Rejected rows remain in their original catalog/FITS index position and are
marked inactive. Candidate totals, size statistics, pairwise chi values,
minimum-chi filtering, connected groups, outputs, and fitting all ignore those
inactive rows. Existing empirical selection thresholds and fitting methods are
unchanged.

Stage 8 no longer repairs NaN ellipticity tokens to `-99`. The Pipeline
Standard reader now rejects non-finite diagnostic tokens; the other live
variants already had no such recovery path.

## Build and test

### Portable Linux/HPC commands

Load a site MPI C++17 compiler plus CFITSIO, FFTW3, Eigen, BLAS, and LAPACK,
then provide their installation prefixes through the existing Make variables:

```sh
make clean
make all CXX=mpicxx STACK_PREFIX="${STACK_PREFIX}" EIGEN_INCLUDE="${EIGEN_INCLUDE}"
make test-psf-candidate-quality CXX=mpicxx STACK_PREFIX="${STACK_PREFIX}" EIGEN_INCLUDE="${EIGEN_INCLUDE}"
./Fourier_Quad_Pipe --help
```

For production, launch the normal pipeline command with the site's MPI or
Slurm launcher. No Windows wrapper or workstation path is required.

### Local validation environment

- Compiler: GCC C++ 15.2.0 through `mpicxx`
- MPI: Open MPI 5.0.10
- CFITSIO: 4.6.4
- FFTW: 3.3.11
- Eigen: 3.4.0
- BLAS/LAPACK: 3.11.0 OpenBLAS build

The local WSL2 override used the existing pixi toolchain and did not install or
modify packages. Standard and Lite each passed a clean full build, the focused
quality-gate regression, and the production `--help` smoke check.

## Validation boundary

The focused regression covers accepted spectra, the zero-median boundary,
NaN and both infinity signs, negative central median with positive total power,
zero/negative total power, invalid stamp shape, and non-finite derived
diagnostics. No real-data Stage 4 -> Stage 5 -> Stage 8 run was available, so
the candidate rejection rate and science-level PSF output still require a
representative dataset check.
