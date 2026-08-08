# process_fd - Field-Distortion Shear Test

## Overview

The fifth pipeline stage (`process_fd`) measures mean galaxy shear as a
function of field distortion. It reads per-exposure shear catalogs
(`*_all.cat`), removes point sources via star-bar fitting, bins sources
by field distortion, and recovers the mean shear per bin with its
uncertainty.

## Feature Switches (FDConfig.hpp)

| Switch | Default | Description |
|--------|---------|-------------|
| `FD_STATIC_MODE` | `PDF_SIGMA` | Statistical mode (see below) |
| `FD_PER_EXPOSURE_STAR_BAR` | `false` | Star-bar mode: `true` = per-exposure fitting (NtoN/HSC), `false` = single global bar (Nto1/DES) |

### `FD_STATIC_MODE` options (enum `StaticMode`)

| Mode | c_best (mean) | sigma (uncertainty) | statis algorithm | jackknife? |
|------|---------------|---------------------|------------------|------------|
| `PDF_SIGMA` | chi2 sign test | quadratic fitting (1/√(2a₁)) | PDF | No |
| `PDF_JACK` | chi2 sign test | jackknife variance | PDF | Yes |
| `SWSE_JACK` | ratio 2·Σ(y·ww)/Σ(ww) | jackknife variance | SWSE | Yes |

## Parameters (FDConfig.hpp = para.inc equivalent)

All parameters are fixed to one set (DES defaults). Key parameters:
- `fd_num=21` spatial bins, `PDF_BINS=4` inner bins, `gf_lim=0.0015`
- `N_jack=50` jackknife regions, `Km_iter=100` k-means iterations
- Quality cuts: `snrlow=20`, `starcut=20`, `chi2_thresh=0.01`, etc.
- Column indices: DES format, 0-based, 44 total columns

## Source Files

| File | Fortran equivalent |
|------|-------------------|
| `process_fd.cpp` | `shear_field_distortion_test_MPI` (main entry) |
| `ShearCatalogReader.cpp` | `read_shear_cat_v2` |
| `StarCutCalculator.cpp` | `calculate_global_star_cut` / `calculate_global_star_cut_auto` / `apply_advanced_cuts` |
| `KMeansClusterer.cpp` | `simple_kmeans_MPI` / `simple_kmeans` / `random_select` |
| `FDMeasurement.cpp` | `plot_comparison_MPI` / `statis_MPI` / `chi2_MPI` / `source_accumulate` |
| `QuadraticFitting.cpp` | `simple_quadratic_fitting` / `determinant3` |

## Usage

```bash
# Run FD test after process_main
mpirun -n 200 Fourier_Quad_Pipe --run-fd true --output-root /path/to/output
```
