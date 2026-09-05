# cpp_Standard

Full C++17 Fourier_Quad pipeline. Choose this variant when a run needs any
optional flat, mask, identity-astrometry, external-PSF, hybrid-PSF, or
PCA/multi-scale branch.

## Build and run

```bash
make clean
make -j4 \
  CXX=/home/alatrion/.pixi/bin/mpicxx \
  STACK_PREFIX=/home/alatrion/.pixi/envs/base \
  EIGEN_INCLUDE=/usr/include/eigen3
cp pipeline.example.ini pipeline.ini
# Edit pipeline.ini with paths and datasets.
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

Standard defaults to `process_init`, `process_main`, `process_rearr`, and
`process_fd` enabled. Review `[process]` before running. `[lensing]` exposes the
Standard branch choices; fixed numerical thresholds remain in
`config/LensingConfig.hpp` and require rebuilding.
Compiled path defaults and output layout names are centralized in
`config/pathconfig.hpp`; INI and CLI values continue to override their runtime
copies.

Archive-format and detector naming conventions are compiled in
`config/InitConfig.hpp`: `ARCHIVE_SUFFIX` selects initializer inputs,
`CCDNUM_KEYWORD` names the physical chip-number FITS keyword for both Science
and DQ HDUs, and
`DQ_STEM_REPLACE_FROM`/`DQ_STEM_REPLACE_TO` map DQ archive stems to science
exposure stems. Their defaults remain `.fits.fz`, `CCDNUM`, and `ood` to `ooi`;
changing any of them requires `make clean` and a rebuild. They intentionally
have no INI or CLI override.

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

`[lensing].astrometry_cat_type` selects legacy large Gaia tiles (`1`) or
1-degree Gaia tiles (`2`); both layouts remain rooted at
`[lensing].astrometry_cat` and use the same RA/Dec row format.

The optional one-time `process_astrocat` phase runs before `process_extcat` and
publishes deduplicated one-degree Gaia tiles. `[astrocat].output_directory` is
independent of `[lensing].astrometry_cat`; configure the consumer path
separately and set `[lensing].astrometry_cat_type = 2` when consuming those
tiles. `PathConfig::ASTROMETRY_TILE_PREFIX` defaults to `astra_`, while
`PathConfig::SOURCE_CAT_TILE_PREFIX` defaults to `extern_` and is shared by
`process_extcat` and external-catalog lookup.

The current tree keeps shared infrastructure in `include/general/` and
`src/general/`, and stage modules under `include/process_*` and
`src/process_*`. Stage 7 writes 28 fields; Stage 9 inserts `EXPO_NUM` immediately
before `ccD_NUM` and appends exposure chi-square. The default external-catalog
row now has 49 fields; regenerate older 48-field products before rearrangement
or FD. Invalid numerical rows are represented by a full `-99999` sentinel row
before Stage 9 rejects them.

Every downstream Norm gate treats an invalid sentinel as an ordinary chip
skip, but a missing or unreadable Norm FITS product is a pipeline-integrity
failure and aborts the MPI world. Stages 4--6 likewise require their Stage-3
catalog files and readable headers; header-only catalogs remain the valid
zero-row representation. Standard PCA residual aggregation schedules physical
CCDNUM directly and opens `<exposure>_<CCDNUM>` StarXY/residual products; it
does not reread Science FITS headers, build or broadcast a positional mapping,
or use a Norm file as an unrelated Stage-6 gate. Missing per-exposure residual
products remain optional aggregation inputs. For external
catalogs, Stage 9 counts all
physical shear/orig lines before production reads, retries one mismatch with
fresh streams, and then consumes exactly the matched data-row count. A
header-only shear catalog remains the zero-source sentinel and is skipped
before the orig file is touched; a zero-line shear file is fatal. Each fixed
iteration reads both paired lines before checking for empty orig content,
missing shear fields, or scientific rejection, and performs no trailing EOF
probe. Shear values retain the fast `stringstream >> float` path: upstream
Stage 7 guarantees no NaN/Inf tokens, so Stage 9 checks only that every row
supplies all `shear_cat_ncols` fields. The Standard non-external path remains
EOF-driven.
`InitializerCCDNUMTest`, `InitializerPublicationCCDNUMTest`, `CCDNUMIdentityTest`,
`PSFExposureTableCCDNUMTest`, and `PSFReconsDirectPathTest` cover extraction
identity failures/resume repair, filename parsing, keyed Astrometry and order
invariance, exposure tables, and direct PCA paths without Science/Norm lookup.
`RequiredInputFailureTest`, `CatalogCombinerLifecycleTest`, and
`PSFReconsOrientationTest` retain adjacent failure/lifecycle/layout coverage.

After the production build, the CCDNUM-focused tests can be reproduced with:

```bash
test_objects=$(find src -name '*.o' \
  ! -path 'src/process_main/MPIScheduler.o' \
  ! -path 'src/process_main/NumericalRecipes.o')
test_cxx=/home/alatrion/.pixi/bin/mpicxx
test_flags='-O2 -std=c++17 -Iinclude -Iconfig -I/usr/include/eigen3 -I/home/alatrion/.pixi/envs/base/include'
test_libs='-L/home/alatrion/.pixi/envs/base/lib -Wl,-rpath,/home/alatrion/.pixi/envs/base/lib -lcfitsio -lfftw3 -lfftw3f -llapack -lblas -lm'

/home/alatrion/.pixi/bin/g++ -O2 -std=c++17 -Iinclude -Iconfig \
  -I/home/alatrion/.pixi/envs/base/include \
  tests/InitializerCCDNUMTest.cpp src/process_init/FitsExtractor.cpp \
  -L/home/alatrion/.pixi/envs/base/lib \
  -Wl,-rpath,/home/alatrion/.pixi/envs/base/lib -lcfitsio \
  -o /tmp/InitializerCCDNUMTest
/tmp/InitializerCCDNUMTest

for test_name in InitializerPublicationCCDNUMTest CCDNUMIdentityTest \
                 PSFExposureTableCCDNUMTest PSFReconsDirectPathTest; do
  test_define=''
  if [ "$test_name" = CCDNUMIdentityTest ]; then
    test_define='-DFQ_STANDARD_VARIANT'
  fi
  $test_cxx $test_flags $test_define "tests/${test_name}.cpp" \
    $test_objects $test_libs -o "/tmp/${test_name}"
  "/tmp/${test_name}"
done
```

See the [main guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
