# cpp_Lite

Reduced C++17 Fourier_Quad pipeline for the fixed production path: Gaia
astrometry, no super-flat, per-chip DQ mask, external source catalog, enabled
deblending, local-polynomial PSF, and no PCA reconstruction. The alternate
Standard branches are absent from this source tree.

## Build and run

```bash
# Edit Initialize.hpp for ordinary compiled defaults.
make clean
make -j4 \
  CXX=/home/alatrion/.pixi/bin/mpicxx \
  STACK_PREFIX=/home/alatrion/.pixi/envs/base \
  EIGEN_INCLUDE=/usr/include/eigen3
cp pipeline.example.ini pipeline.ini
# Optionally override defaults for this run in pipeline.ini.
/home/alatrion/.pixi/bin/mpirun -np 4 \
  ./Fourier_Quad_Pipe --config pipeline.ini
```

For a nonstandard library prefix, pass `STACK_PREFIX`; pass `EIGEN_INCLUDE`
when Eigen is elsewhere. The current Makefile exposes only the `all` and
`clean` targets. Focused tests are compiled explicitly rather than through
named Make targets.

The verified WSL toolchain is GCC 15.2.0 through Open MPI 5.0.10, CFITSIO 4.6.4,
FFTW 3.3.11, Eigen 3.4.0, and OpenBLAS-backed BLAS/LAPACK 3.11.0. The commands
above select the pixi compiler/library prefix explicitly; other sites must
provide equivalent C++17 MPI, Eigen3, CFITSIO, FFTW3, LAPACK, and BLAS
dependencies.

Lite defaults to `process_init` and `process_main` enabled, with
`process_rearr` and `process_fd` disabled. It rejects Standard-only `[lensing]`
keys rather than silently ignoring them.
Ordinary compiled defaults are centralized in `Initialize.hpp`; `config/`
retains advanced/internal defaults, and INI and CLI values continue to override
their runtime copies. Lite always consumes repartitioned one-degree Gaia tiles
rooted at `[lensing].astrometry_cat`.

## CCDNUM-native output contract

Initializer manifest schema 3 names every extracted chip
`<exposure>_<CCDNUM>.fits` for both Science and DQ. A missing `CCDNUM` skips
only that two-dimensional HDU; malformed, non-positive, or duplicate values
fail the archive. Resume accepts an existing chip only when its FITS header
matches the CCDNUM encoded in the filename. Exposure-list row order is not a
chip identity and may change after partial resume.

All downstream persistent interfaces use physical CCDNUM: Astrometry `.head`
rows are keyed lookups; StarInfo, StarComp, ExpoInfo, Stage-9 `ccD_NUM`, and
chip-local paths keep that key. Common Astrometry CRVAL is chosen from the
smallest initially valid CCDNUM, so list reordering cannot change the
projection center. This is an intentional schema break: start with a clean
output root and rerun `process_init`; old continuously numbered trees are not
supported or auto-detected.

The optional one-time `process_astrocat` phase runs before `process_extcat` and
publishes deduplicated one-degree Gaia tiles. `[astrocat].output_directory` is
independent of `[lensing].astrometry_cat`; configure producer and consumer paths
separately. Likewise, `[extcat].output_directory` is the producer destination,
while `[lensing].source_cat` is the `process_main` consumer input. The compiled
tile prefixes default to `astra_` and `extern_`.

The current tree keeps shared infrastructure in `include/general/` and
`src/general/`. Stage 7 writes 28 fields; Stage 9 inserts `EXPO_NUM` immediately
before `ccD_NUM` and appends exposure chi-square. The default Lite row now has
49 fields; regenerate older 48-field products before rearrangement or FD.
Full-row `-99999` sentinels remain the invalid numerical-source contract.

Every downstream Norm gate treats an invalid sentinel as an ordinary chip
skip, but a missing or unreadable Norm FITS product is a pipeline-integrity
failure and aborts the MPI world. Stages 4--6 likewise require their Stage-3
catalog files and readable headers; header-only catalogs remain the valid
zero-row representation. Stage 9 counts all physical shear/orig lines
before production reads, retries one mismatch with fresh streams, and then
consumes exactly the matched data-row count. A header-only shear catalog
remains the zero-source sentinel and is skipped before the orig file is
touched; a zero-line shear file is fatal. Each fixed iteration reads both
paired lines before checking for empty orig content, missing shear fields, or
scientific rejection, and performs no trailing EOF probe. Shear values retain
the fast `stringstream >> float` path: upstream Stage 7 guarantees no NaN/Inf
tokens, so Stage 9 checks only that every row supplies all
`shear_cat_ncols` fields.
`InitializerCCDNUMTest`, `InitializerPublicationCCDNUMTest`,
`CCDNUMIdentityTest`, and
`PSFExposureTableCCDNUMTest` cover extraction identity failures/resume repair,
filename parsing, keyed Astrometry and order invariance, and Stage-5 exposure
tables. `RequiredInputFailureTest` and `CatalogCombinerLifecycleTest` retain
adjacent failure/lifecycle coverage.

After the production build, the CCDNUM-focused tests can be reproduced with:

```bash
test_objects=$(find src -name '*.o' \
  ! -path 'src/process_main/MPIScheduler.o' \
  ! -path 'src/process_main/NumericalRecipes.o')
test_cxx=/home/alatrion/.pixi/bin/mpicxx
test_flags='-O2 -std=c++17 -I. -Iinclude -Iconfig -I/usr/include/eigen3 -I/home/alatrion/.pixi/envs/base/include'
test_libs='-L/home/alatrion/.pixi/envs/base/lib -Wl,-rpath,/home/alatrion/.pixi/envs/base/lib -lcfitsio -lfftw3 -lfftw3f -llapack -lblas -lm'

/home/alatrion/.pixi/bin/g++ -O2 -std=c++17 -I. -Iinclude -Iconfig \
  -I/home/alatrion/.pixi/envs/base/include \
  tests/InitializerCCDNUMTest.cpp src/process_init/FitsExtractor.cpp \
  -L/home/alatrion/.pixi/envs/base/lib \
  -Wl,-rpath,/home/alatrion/.pixi/envs/base/lib -lcfitsio \
  -o /tmp/InitializerCCDNUMTest-Lite
/tmp/InitializerCCDNUMTest-Lite

for test_name in InitializerPublicationCCDNUMTest CCDNUMIdentityTest \
                 PSFExposureTableCCDNUMTest; do
  $test_cxx $test_flags "tests/${test_name}.cpp" \
    $test_objects $test_libs -o "/tmp/${test_name}-Lite"
  "/tmp/${test_name}-Lite"
done
```

See the [main guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
