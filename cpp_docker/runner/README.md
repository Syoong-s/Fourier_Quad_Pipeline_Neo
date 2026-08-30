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

Set the SIF or OCI/archive source, `CPP_SOURCE_HOST`, `PROCESS_DATA_HOST`,
catalog/calibration binds, and Apptainer cache/tmp paths. Set Science, DQ,
extcat input, rearr output, exposure-list, and FD output binds only when the
selected phases require them.

`cpppipeline.env` is sourced as Bash. Keep `HPC_MODULES`,
`HPC_EXTRA_BINDS`, `HPC_PASSTHROUGH_ENV`, `HPC_CONTAINER_ENV`, and
`SRUN_ARGS` as indexed arrays. Preserve:

```text
MPI_LAUNCH_MODE=srun
SLURM_MPI_TYPE=pmi2
```

Modules may expose Slurm or Apptainer, but must not inject host MPI libraries
into the application.

## Acquire the SIF

From a reviewed Docker archive:

```bash
sbatch build-sif.slurm
```

From a digest-pinned registry image:

```bash
bash pull-sif.sh
```

Both paths refuse to overwrite an existing SIF and create a SHA256 sidecar.
Set `CPP_SIF_SHA256_EXPECTED` from that sidecar before production.

## Validate and run

```bash
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
sbatch mpi-smoke-test.slurm
```

Only after these checks, launch the pipeline:

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
