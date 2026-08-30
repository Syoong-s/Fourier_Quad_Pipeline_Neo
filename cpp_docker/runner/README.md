# Fourier_Quad C++ on Slurm/Apptainer

This runner launches the container-linked pipeline on x86_64 Slurm clusters
that advertise `pmi2`. The SIF supplies the compiler, OpenMPI, and scientific
libraries; host compiler/OpenMPI ABI compatibility is not required.

> 中文版：[README-CN.md](README-CN.md)

## Prerequisites

- `pmi2` appears in `srun --mpi=list`;
- Apptainer or Singularity is available on compute nodes;
- runner, SIF, source, catalogs, data, cache, and temporary directories are on
  storage visible at the same paths from every node;
- allocated nodes can communicate over TCP.

Run the read-only site check first:

```bash
bash inspect-cluster-mpi.sh
```

## Configure

```bash
cp cpppipeline.env.example cpppipeline.env
```

### Common runner parameters

Host paths in `cpppipeline.env` must be visible at the same location on every allocated
node. Paired container paths must match the paths used by `pipeline.ini` or the
corresponding CLI options.

| Parameter or location | Typical change | Constraint |
|---|---|---|
| `OCI_IMAGE_URI` / `CPP_DOCKER_ARCHIVE` | Choose either OCI pull or a local Docker archive build. | Pin `OCI_IMAGE_URI` by digest; an archive must be an x86_64 image visible to the compute environment. |
| `CPP_SIF`, `CPP_SIF_SHA256_EXPECTED` | Set the shared SIF path and optional expected SHA256. | Keep one path across all nodes; fill the checksum after the sidecar is generated. |
| `CPP_SOURCE_HOST/CONTAINER` | Point to `cpp_Standard` or `cpp_Lite` and `/workspace/src_pipe` in the container. | Host source must be writable for `pipeline.ini` and build products. |
| `SCIENCE_ROOT_*`, `DQ_ROOT_*` | Set the Science/DQ archives used by initialization/main. | Bind only when needed; container paths must match `[init]` or CLI, and Lite main requires per-CCD DQ. |
| `ASTROMETRY_CAT_*`, `SOURCE_CAT_*`, `FLAT_PATH_*` | Point to Gaia, external-source, and flat directories. | Container paths must match effective `[lensing]`/`[extcat]` paths; flat is needed only for its Standard branch. |
| `PROCESS_DATA_*`, `CPP_EXPO_LIST_CONTAINER` | Set shared writable processing storage and the default exposure-list container path. | Exposure lists and configuration must lie below a bound container path. |
| `EXTCAT_INPUT_*`, `REARR_OUTPUT_*`, `EXPOLIST_DIR_*`, `FD_OUTPUT_*` | Set only for phases that need independent directories. | Unset fields are not bound; INI/CLI must use the corresponding container paths. |
| `HPC_SHARED_SCRATCH_HOST`, `APPTAINER_CACHE_DIR`, `APPTAINER_TMP_DIR` | Set site shared writable scratch/cache/tmp. | Provide enough space and compute-node visibility. |
| `APPTAINER_BIN`, `HPC_MODULES`, `SITE_ENV_SCRIPT` | Adapt to site commands and module setup. | Keep `HPC_MODULES` as a Bash array and do not inject host MPI libraries. |
| `MPI_LAUNCH_MODE=srun`, `SLURM_MPI_TYPE=pmi2`, `SRUN_ARGS=()` | Preserve the launcher and add site flags to `SRUN_ARGS` only when required. | This runner requires `srun` + `pmi2`; do not flatten the array to a scalar string. |
| `HPC_EXTRA_BINDS=()`, `HPC_PASSTHROUGH_ENV`, `HPC_CONTAINER_ENV=()` | Set only for genuinely required extra directories or environment values. | Keep all three as Bash arrays; the default passthrough contains only `OMP_NUM_THREADS`, and host state passed through `--cleanenv` should remain minimal. |
| `CPP_BUILD_JOBS`, `CPP_MAKE_CLEAN`, `CPP_EXECUTABLE` | Adjust build parallelism, clean-before-build, and executable container path. | Never let concurrent jobs clean or compile the same source copy. |
| `CPP_IMAGE_ID_EXPECTED` and compiler/MPI expected versions | Update only with a changed image stack. | Keep the expectations consistent with actual image versions. |
| `#SBATCH` directives in `cpppipeline.slurm` | Set partition, nodes, tasks, tasks-per-node, CPU, memory, time, and log paths. | Follow site policy and keep task counts consistent with the MPI-rank plan. |

`cpppipeline.env` is sourced as Bash. Keep `HPC_MODULES`,
`HPC_EXTRA_BINDS`, `HPC_PASSTHROUGH_ENV`, `HPC_CONTAINER_ENV`, and
`SRUN_ARGS` as indexed arrays. Preserve:

```text
MPI_LAUNCH_MODE=srun
SLURM_MPI_TYPE=pmi2
```

Modules may expose Slurm or Apptainer, but must not inject host MPI libraries
into the application.

## Acquire and validate

Use one acquisition path:

```bash
sbatch build-sif.slurm          # existing local Docker archive
bash pull-sif.sh                # download an OCI image from a remote registry
```

Both paths refuse to overwrite an existing SIF and create a SHA256 sidecar.
Set `CPP_SIF_SHA256_EXPECTED` from that sidecar before production, then we
recommend running:

```bash
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
sbatch mpi-smoke-test.slurm
```

## Run

Pass an INI path or C++ CLI options after the script name. For example:

```bash
sbatch cpppipeline.slurm \
  --config /workspace/src_pipe/pipeline.ini
```

Arguments after the script name are passed unchanged to
`Fourier_Quad_Pipe`. With no arguments, the runner passes
`CPP_EXPO_LIST_CONTAINER` as the legacy exposure-list argument.

The source bind is writable because compilation products stay there. Do not
compile the same source copy concurrently. Catalog/calibration binds are
read-only and processing/output binds are writable.

Slurm resource directives are templates. Override partition, account, nodes,
tasks, CPUs, memory, time, and logs according to the target site without
changing the `srun --mpi=pmi2` process boundary.

For `process_astrocat`, expose the raw Gaia directory as a read-only entry in
`HPC_EXTRA_BINDS`, then write the generated tiles below the writable
`PROCESS_DATA_CONTAINER` or another explicit writable bind. Use those
container paths in `[astrocat]` or the corresponding CLI options. The producer
output is independent of `[lensing].astrometry_cat`; configure that consumer
path and `[lensing].astrometry_cat_type = 2` separately.
