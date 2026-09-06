# cpp_Lite

Reduced C++17 Fourier_Quad pipeline for the fixed production path: Gaia
astrometry, no super-flat, per-chip DQ masks, external source catalogs,
deblending, local-polynomial PSF, and no PCA reconstruction. Standard-only
selector branches are absent from this source tree.

## Configure, build, and run

Set the ordinary compiled defaults in `Initialize.hpp` first:

```bash
make -j4
cp pipeline.example.ini pipeline.ini
# Optionally override defaults for this run in pipeline.ini.
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe --config pipeline.ini
```

Pass `CXX`, `STACK_PREFIX`, and `EIGEN_INCLUDE` when the MPI compiler,
scientific libraries, or Eigen headers are outside default search paths. The
current Makefile exposes only `all` and `clean`; focused test sources under
`tests/` are compiled explicitly when needed.

The required stack is an MPI C++ compiler with C++17 support, CFITSIO, FFTW3
(double and single precision), Eigen3, LAPACK, and BLAS. The published
container pins GCC 12.3.0, OpenMPI 4.1.8, CFITSIO 4.6.4, FFTW 3.3.11, Eigen
3.4.0, LAPACK 3.11.0, and OpenBLAS 0.3.33.

## Configuration

Lite defaults to initialization and main processing enabled, with rearrangement
and FD disabled. `Initialize.hpp` centralizes ordinary instrument, workflow,
dataset, path, catalog-schema, archive, and FITS naming defaults.
`pipeline.example.ini` provides per-run overrides, with command-line options
taking final precedence:

```text
Initialize.hpp and other compiled defaults < --config INI < command-line options
```

Lite rejects Standard-only `[lensing]` keys rather than silently ignoring
them. It always uses per-chip DQ masks and repartitioned one-degree Gaia tiles.
`process_astrocat` and `process_extcat` producer outputs must be configured
separately from the catalog paths consumed by `process_main`.

## Data compatibility

Initializer manifest schema 3 names Science and DQ products with physical
`CCDNUM`, and downstream interfaces use the same identity rather than
exposure-list position. This is incompatible with older continuously numbered
output trees; start with a clean output root and rerun initialization.

Stage 7 writes 28 fields and Stage 9 adds exposure chi-square. With the default
18-field external catalog and two identity fields, the final row has 49 fields;
regenerate older 48-field products before rearrangement or FD. Header-only
catalogs remain the valid zero-row representation, while missing or unreadable
required products stop the MPI job.

See the [C++ guide](../CPP_GUIDE.md) for phase, I/O, schema, and failure
contracts, and the
[parameter reference](../CPP_PIPELINE_PARAMETERS.md) for all INI, CLI, and
compile-time settings.
