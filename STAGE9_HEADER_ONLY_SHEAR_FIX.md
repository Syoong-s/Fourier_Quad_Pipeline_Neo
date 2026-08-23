# Stage-9 header-only shear handling

## Behavior

Stage 7 uses a valid catalog header with zero data rows as the chip-level
sentinel for an invalid PSF, failed astrometry, or another no-shear condition.
Stage 9 now classifies and skips that sentinel before resolving or counting the
paired `_orig.cat` file. A nonempty shear catalog still invokes
`requireMatchingCatalogDataRows`, so every physical-row mismatch remains fatal.

The correction is implemented independently in `cpp_Standard` and `cpp_Lite`.
It does not change catalog columns, calibration, source cuts, PSF behavior, or
the strict row-count helper.

## Regression coverage

Each variant has `tests/CatalogCombinerLifecycleTest.cpp`. Its sentinel case
writes a header-only shear catalog together with a populated original catalog
and requires Stage 9 to remove stale output without creating a replacement.
The existing `CatalogRowCountTest.cpp` continues to require unequal nonempty
catalogs to terminate through the production fail-fast path.

Because this worktree currently has local Makefile changes that omit standalone
test targets, the lifecycle and row-count tests were compiled directly against
the production sources. A portable lifecycle-test command from either variant
directory is:

```bash
mpicxx -Iinclude -Iconfig -Iinclude/process_extcat -Iinclude/process_init \
  -Iinclude/process_main -Iinclude/process_rearr -Iinclude/process_fd \
  -std=c++17 tests/CatalogCombinerLifecycleTest.cpp \
  src/process_main/CatalogCombiner.cpp src/process_main/Universalblock.cpp \
  src/process_main/UniversalUtils.cpp src/process_main/FitsIO.cpp \
  src/process_main/LinearSolve.cpp -o test_catalog_lifecycle \
  -lcfitsio -lfftw3 -lfftw3f -llapack -lblas -lm
```

For Standard, insert `src/RuntimeConfig.cpp` before `-o`; Lite does not require
that additional translation unit.

## Build and run environments

- Local verification: GCC/mpicxx 15.2.0, Open MPI 5.0.10, CFITSIO 4.6.3,
  FFTW 3.3.10, Eigen 3.4.0, LAPACK/BLAS from the existing Pixi science stack.
- Local override:

  ```bash
  make -j2 CXX=/home/alatrion/.pixi/bin/mpicxx \
    STACK_PREFIX=/home/alatrion/.pixi/envs/base \
    EIGEN_INCLUDE=/usr/include/eigen3
  ```

- Portable cluster build: load the site GCC, MPI, CFITSIO, FFTW, Eigen, and
  LAPACK/BLAS modules, then run `make -j2 CXX=mpicxx` with `STACK_PREFIX` and
  `EIGEN_INCLUDE` set only when those modules are outside compiler defaults.
- Local smoke: `./Fourier_Quad_Pipe --help`.
- Cluster continuation: rebuild the selected variant and run Stage 8 plus
  Stage 9 with `PROCESS_stage=19*23`; the current driver rejects Stage 9 alone.
