# cpp_Standard

Full C++17 Fourier_Quad pipeline. Choose this variant when a run needs any
optional flat, mask, identity-astrometry, external-PSF, hybrid-PSF, or
PCA/multi-scale branch.

## Build and run

```bash
make -j4
cp pipeline.example.ini pipeline.ini
# Edit pipeline.ini with paths and datasets.
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

Standard defaults to `process_init`, `process_main`, `process_rearr`, and
`process_fd` enabled. Review `[process]` before running. `[lensing]` exposes the
Standard branch choices; fixed numerical thresholds remain in
`config/LensingConfig.hpp` and require rebuilding.
Compiled path defaults and output layout names are centralized in
`config/pathconfig.hpp`; INI and CLI values continue to override their runtime
copies.

Archive-format and detector naming conventions are compiled in
`config/InitConfig.hpp`: `ARCHIVE_SUFFIX` selects initializer inputs,
`CCDNUM_KEYWORD` names the DQ/main chip-number FITS keyword, and
`DQ_STEM_REPLACE_FROM`/`DQ_STEM_REPLACE_TO` map DQ archive stems to science
exposure stems. Their defaults remain `.fits.fz`, `CCDNUM`, and `ood` to `ooi`;
changing any of them requires `make clean` and a rebuild. They intentionally
have no INI or CLI override.

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
failure and aborts the MPI world. For external catalogs, Stage 9 counts all
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
`UniversalblockTest` and `CatalogCombinerLifecycleTest` provide focused local
coverage for these contracts and are compiled explicitly with the same C++17
MPI and science-library settings as the production build.

See the [main guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
