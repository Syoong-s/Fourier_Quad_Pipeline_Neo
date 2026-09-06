# Fourier_Quad C++ on Slurm/Apptainer

This runner launches the container-linked Pipeline on x86_64 Slurm clusters
that expose `pmi2`. The SIF supplies the compiler, OpenMPI, and scientific
libraries; host compiler/OpenMPI ABI compatibility is not required.

> 中文版：[README-CN.md](README-CN.md)

## Prerequisites and configuration

Confirm that `pmi2` appears in `srun --mpi=list`, Apptainer/Singularity
works on compute nodes, allocated nodes can communicate over TCP, and all
paths are visible at identical locations from every node. Then run:

```bash
bash inspect-cluster-mpi.sh
cp cpppipeline.env.example cpppipeline.env
```

### Common runner parameters

Host paths in `cpppipeline.env` must be shared across allocated nodes. Paired
container paths must match `pipeline.ini` or the corresponding CLI options.

| Parameter or location | Typical setting | Constraint |
|---|---|---|
| `OCI_IMAGE_URI` / `CPP_DOCKER_ARCHIVE` | Choose OCI pulling or conversion of a local Docker archive. | Pin the URI by digest; an archive must be x86_64 and visible to the build job. |
| `CPP_SIF`, `CPP_SIF_SHA256_EXPECTED` | Shared SIF path and optional expected SHA256. | Use one path on all nodes and fill the checksum after acquisition. |
| `CPP_SOURCE_HOST/CONTAINER` | `cpp_Standard` or `cpp_Lite`, normally mounted at `/workspace/src_pipe`. | Host source must be writable so configuration and build products persist. |
| `SCIENCE_ROOT_*`, `DQ_ROOT_*` | Science/DQ archives used by initialization or main processing. | Bind only when needed; Lite main processing requires per-chip DQ. |
| `ASTROMETRY_CAT_*`, `SOURCE_CAT_*`, `FLAT_PATH_*` | Gaia, external-source, and flat directories. | Container paths must match effective INI/CLI values; flat is Standard-only. |
| `PROCESS_DATA_*`, `CPP_EXPO_LIST_CONTAINER` | Shared writable processing directory and default exposure-list path. | Exposure lists and configuration must resolve below a bind. |
| `EXTCAT_INPUT_*`, `REARR_OUTPUT_*`, `EXPOLIST_DIR_*`, `FD_OUTPUT_*` | Optional independent phase directories. | Unset pairs are not bound; INI/CLI must use the matching container paths. |
| `HPC_SHARED_SCRATCH_HOST`, `APPTAINER_CACHE_DIR`, `APPTAINER_TMP_DIR` | Shared writable scratch/cache/tmp. | Must have sufficient space and be compute-node accessible. |
| `APPTAINER_BIN`, `HPC_MODULES`, `SITE_ENV_SCRIPT` | Site runtime and module setup. | Do not inject host MPI libraries. |
| `MPI_LAUNCH_MODE=srun`, `SLURM_MPI_TYPE=pmi2`, `SRUN_ARGS=()` | Preserve the launch mode; add site flags through `SRUN_ARGS` if needed. | This runner requires `srun` + `pmi2`; keep array syntax. |
| `HPC_EXTRA_BINDS=()`, `HPC_PASSTHROUGH_ENV=(OMP_NUM_THREADS)`, `HPC_CONTAINER_ENV=()` | Required extra files or environment values only. | Keep Bash arrays and minimize state passed through `--cleanenv`. |
| `CPP_BUILD_JOBS`, `CPP_MAKE_CLEAN`, `CPP_EXECUTABLE` | Build parallelism, clean-before-build, and executable path. | Do not let concurrent jobs clean or compile the same source copy. |
| Image/version expectations | Update with a changed image stack. | Must match the versions actually present in the SIF. |
| `#SBATCH` in `cpppipeline.slurm` | Site partition, account, nodes, tasks, CPU, memory, time, and logs. | Follow site policy and keep task counts consistent with MPI ranks. |

`cpppipeline.env` is Bash. Preserve its indexed arrays and keep:

```text
MPI_LAUNCH_MODE=srun
SLURM_MPI_TYPE=pmi2
```

Modules may expose Slurm or Apptainer but must not inject host MPI libraries
into the application.

## Acquire and validate

Use one acquisition path:

```bash
sbatch build-sif.slurm          # existing local Docker archive
bash pull-sif.sh                # remote OCI image
```

Both paths refuse to overwrite an existing SIF and create a SHA256 sidecar.
Set `CPP_SIF_SHA256_EXPECTED` from that sidecar before production, then run:

```bash
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
sbatch mpi-smoke-test.slurm
```

## Run

Pass an INI path or C++ CLI options after the script name:

```bash
sbatch cpppipeline.slurm \
  --config /workspace/src_pipe/pipeline.ini
```

Arguments are forwarded unchanged to `Fourier_Quad_Pipe`. With no arguments,
the runner passes `CPP_EXPO_LIST_CONTAINER` as the legacy exposure-list
argument. The source bind is writable for compilation; do not build the same
copy concurrently. Adjust Slurm resources for the site without changing the
`srun --mpi=pmi2` launch boundary.

For `process_astrocat`, add the raw Gaia directory to `HPC_EXTRA_BINDS` as
a read-only bind and write generated tiles below a writable bind. Configure
the producer paths in `[astrocat]` or CLI; configure the independent
`[lensing].astrometry_cat` consumer path separately. Standard must also set
`[lensing].astrometry_cat_type = 2`; Lite already uses the one-degree layout.
