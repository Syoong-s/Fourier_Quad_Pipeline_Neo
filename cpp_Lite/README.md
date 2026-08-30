# cpp_Lite

Reduced C++17 Fourier_Quad pipeline for the fixed production path: Gaia
astrometry, no super-flat, per-chip DQ mask, external source catalog, enabled
deblending, local-polynomial PSF, and no PCA reconstruction. The alternate
Standard branches are absent from this source tree.

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

Lite defaults to `process_init` and `process_main` enabled, with
`process_rearr` and `process_fd` disabled. It rejects Standard-only `[lensing]`
keys rather than silently ignoring them.
`[lensing].astrometry_cat_type` selects legacy large Gaia tiles (`1`) or
1-degree Gaia tiles (`2`); both layouts remain rooted at
`[lensing].astrometry_cat` and use the same RA/Dec row format.

The optional one-time `process_astrocat` phase runs before `process_extcat` and
publishes deduplicated one-degree Gaia tiles. `[astrocat].output_directory` is
independent of `[lensing].astrometry_cat`; configure the consumer path
separately and set `[lensing].astrometry_cat_type = 2` when consuming those
tiles.

The current tree keeps shared infrastructure in `include/general/` and
`src/general/`. Stage 7 writes 28 fields; Stage 9 inserts `EXPO_NUM` immediately
before `ccD_NUM` and appends exposure chi-square. The default Lite row now has
49 fields; regenerate older 48-field products before rearrangement or FD.
Full-row `-99999` sentinels remain the invalid numerical-source contract.

See the [main guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
