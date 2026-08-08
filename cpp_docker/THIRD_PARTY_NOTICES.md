# Third-party notices

This directory does not redistribute downloaded dependency source archives.
Docker builds install binary packages from conda-forge according to
`pixi.lock`, use digest-pinned Pixi and Rocky Linux images, and build the
Slurm 25.11.2 PMI2 client plus OpenMPI 4.1.8 from the upstream archives
recorded in `SOURCES.md`.

GCC, binutils, OpenMPI, Slurm, CFITSIO, FFTW, Eigen, LAPACK, BLAS, OpenBLAS,
Rocky Linux, Pixi, and all transitive packages remain subject to their
respective upstream licenses. `pixi.lock` records package-level license
metadata where supplied by conda-forge. Review the license files shipped in
the final OCI image or SIF and the upstream license terms before
redistributing a built artifact.

Repository-authored Docker, runner, and documentation files do not replace or
modify third-party license obligations. A SIF produced from the image is a
binary redistribution of those components and must retain the applicable
notices.
