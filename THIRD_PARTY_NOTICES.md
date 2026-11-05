# 📜 Third-Party Notices

[🏠 Home](README.md) · [⚖️ MIT License](LICENSE)

This repository does not redistribute downloaded third-party source archives. Container toolchain images (`cpp_docker/Dockerfile`) install binary packages from `conda-forge` locked via `cpp_docker/pixi.lock`, use digest-pinned `Pixi` and `Rocky Linux` base images, and build the Slurm PMI2 client (`25.11.2`) and OpenMPI (`4.1.8`) from upstream archives verified by SHA-256 checksums.

| Component Category | Covered Upstream Packages | License Obligation |
| :--- | :--- | :--- |
| **Compilers & Build Tools** | GCC (`12.3.0`), binutils, GNU Make, Pixi | Governed by respective upstream open-source licenses (`pixi.lock`). |
| **MPI & Scheduler Clients** | OpenMPI (`4.1.8`), Slurm PMI2 client (`25.11.2`) | Subject to BSD / GPL / upstream terms recorded in `Dockerfile`. |
| **Scientific Libraries** | CFITSIO (`4.6.4`), FFTW (`3.3.11`), Eigen (`3.4.0`), OpenBLAS, LAPACK, BLAS | Retain upstream licenses (MIT, MPL-2.0, GPL, BSD, NASA/CFITSIO). |
| **Base Operating System** | Rocky Linux (`8.10`) and transitive dependencies | Subject to upstream distribution licensing terms. |

> [!IMPORTANT]
> Repository-authored C++ source, Docker tooling, HPC runner scripts, and documentation (covered by the [MIT License](LICENSE)) do not replace or modify third-party license obligations. Any OCI container image or Apptainer `.sif` file built from `cpp_docker/` constitutes a binary redistribution of those upstream components and must retain all applicable notices.
