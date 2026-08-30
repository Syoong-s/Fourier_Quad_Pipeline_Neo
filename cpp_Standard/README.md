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
`[lensing].astrometry_cat_type` selects legacy large Gaia tiles (`1`) or
1-degree Gaia tiles (`2`); both layouts remain rooted at
`[lensing].astrometry_cat` and use the same RA/Dec row format.

The optional one-time `process_astrocat` phase runs before `process_extcat` and
publishes deduplicated one-degree Gaia tiles. `[astrocat].output_directory` is
independent of `[lensing].astrometry_cat`; configure the consumer path
separately and set `[lensing].astrometry_cat_type = 2` when consuming those
tiles.

The current tree keeps shared infrastructure in `include/general/` and
`src/general/`, and stage modules under `include/process_*` and
`src/process_*`. Stage 7 writes 28 fields; Stage 9 inserts `EXPO_NUM` immediately
before `ccD_NUM` and appends exposure chi-square. The default external-catalog
row now has 49 fields; regenerate older 48-field products before rearrangement
or FD. Invalid numerical rows are represented by a full `-99999` sentinel row
before Stage 9 rejects them.

Every downstream Norm gate treats an invalid sentinel as an ordinary chip
skip, but a missing or unreadable Norm FITS product is a pipeline-integrity
failure and aborts the MPI world. Stage 9 likewise aborts if a paired external
catalog row cannot be read after row-count preflight; it never treats runtime
EOF as normal completion. `UniversalblockTest`, `CatalogRowCountTest`, and
`CatalogCombinerLifecycleTest` provide the focused local regression coverage
for these contracts and are compiled explicitly with the same C++17 MPI and
science-library settings as the production build.

See the [main guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
