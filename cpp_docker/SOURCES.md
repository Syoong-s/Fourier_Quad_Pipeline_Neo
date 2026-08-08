# Source provenance and integrity

The Docker build uses two digest-pinned images:

| Purpose | Image |
| --- | --- |
| Pixi executable | `ghcr.io/prefix-dev/pixi:0.71.1@sha256:cc1401e1553164ad032908e1f6a1dbbaa8c8e067ce4347ae831c43a066bbc66c` |
| Runtime and build base | `docker.io/rockylinux/rockylinux:8.10@sha256:f5529992e67440c1a4ae7788244d4381c6909159a88eacd95b7523ae47ced82e` |

The Rocky digest is the verified linux/amd64 manifest. The Dockerfile rejects
non-amd64 builds.

`pixi.lock` contains one environment named `default`. It records conda-forge
package URLs, builds, dependency graphs, archive checksums, and license
metadata. `pixi install --locked` rejects an unresolved or changed
environment. MPI and PMIx are intentionally absent from that environment.

The Dockerfile also downloads and verifies:

| Component | Version | Upstream archive | SHA-256 |
| --- | --- | --- | --- |
| Slurm PMI2 client source | 25.11.2 | <https://download.schedmd.com/slurm/slurm-25.11.2.tar.bz2> | `9bfd844f746c1268d6e0cb14bda586b548b1c4b12f8c72181397d0cec9ce5d39` |
| OpenMPI | 4.1.8 | <https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.8.tar.bz2> | `466f68e3132a1dc02710cc2011fafced8336d98359fa2dae4dddcfd5719f12a9` |

Only Slurm's PMI2 client library and header are installed. The image does not
contain a Slurm controller, daemon, scheduler command set, MUNGE, or a Slurm
configuration. OpenMPI is built against that PMI2 client with the locked
G++ 12.3.0 toolchain. OpenMPI's own Slurm direct-launch environment components
are enabled, but they do not add a dependency on host `libslurm`.

No downloaded source archive or unpacked dependency source tree is copied into
the final runtime stage. Pipeline source and science data are never copied into
the image.
