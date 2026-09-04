# cpp_Lite

Reduced C++17 Fourier_Quad pipeline for the fixed production path: Gaia
astrometry, no super-flat, per-chip DQ mask, external source catalog, enabled
deblending, local-polynomial PSF, and no PCA reconstruction. The alternate
Standard branches are absent from this source tree.

## Build and run

```bash
# Edit Initialize.hpp for ordinary compiled defaults.
make -j4
cp pipeline.example.ini pipeline.ini
# Optionally override defaults for this run in pipeline.ini.
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

For a nonstandard library prefix, pass `STACK_PREFIX`; pass `EIGEN_INCLUDE`
when Eigen is elsewhere. The current Makefile exposes only the `all` and
`clean` targets. Focused tests are compiled explicitly rather than through
named Make targets.

Local focused verification uses the MPI C++ wrapper from GCC 15.2.0, with
CFITSIO 4.6.3 and FFTW3 3.3.10 available. The local full build uses Eigen3 from
`/usr/include/eigen3`; other sites must provide equivalent C++17 MPI, Eigen3,
LAPACK, and BLAS dependencies.

Lite defaults to `process_init` and `process_main` enabled, with
`process_rearr` and `process_fd` disabled. It rejects Standard-only `[lensing]`
keys rather than silently ignoring them.
Ordinary compiled defaults are centralized in `Initialize.hpp`; `config/`
retains advanced/internal defaults, and INI and CLI values continue to override
their runtime copies. Lite always consumes repartitioned one-degree Gaia tiles
rooted at `[lensing].astrometry_cat`.

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
`UniversalblockTest`, `RequiredInputFailureTest`, and
`CatalogCombinerLifecycleTest` provide focused local coverage for these
contracts and are compiled explicitly with the same C++17 MPI and
science-library settings as the production build.

See the [main guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
