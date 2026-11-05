# HPC with Apptainer

[Home](README.md) · [C++ user guide](CPP_GUIDE.md) · [Container paths](PATHS.md)

For native builds, use the [C++ guide](CPP_GUIDE.md) and your site's MPI launcher.
This page covers the optional container runner. It targets x86_64 systems with
Slurm PMI2 and an Apptainer/Singularity runtime supporting `--cleanenv`,
`--no-eval` and `--no-mount`. The supplied toolchain image uses its own OpenMPI
4.1.8 and PMI2 client. Compatibility with the site's Slurm/MPI network is a site
requirement, not something a configuration check can certify.

## Prepare once

Use `cpp_docker/fqctl.py init`, edit the JSON, then run `configure` as described
in [Docker configuration](DOCKER.md#1-configure-host-paths). These steps need
Python and shared filesystems, but **do not require Docker**. Set the following
HPC fields before generating:

| JSON field | Example / meaning |
| --- | --- |
| `hpc.sif` | `/shared/images/cpppipeline.sif`; an image readable by every compute node. |
| `hpc.runtime` | `apptainer`, `singularity`, or its executable path. |
| `hpc.srun_args` | Optional placement settings, e.g. `["--cpu-bind=cores"]`. |
| `hpc.environment` | Additional container values, e.g. `{"OMP_NUM_THREADS":"1"}`. |

Source, JSON-generated configuration, data and the runner itself must be
accessible from the participating nodes. Load site modules in your own shell
or job script. Do not bind a host MPI over the container's MPI.

Obtain a compatible SIF using your site's supported image-transfer procedure.
For example, on a node allowed to download images:

```bash
apptainer pull /shared/images/cpppipeline.sif \
  docker://ghcr.io/syoong-s/fourier_quad_pipeline_neo:gxx12.3-openmpi4.1.8-pmi2
```

The image is separate from the source and scientific data. No image download
or conversion is hidden inside the run command.

## Check configuration

```bash
bash /absolute/repo/cpp_docker/runner/run.sh --check \
  /absolute/run/cpppipeline.env
```

This explicitly runs the configuration check. It does not start a container,
require an allocation or require the SIF to have been downloaded. It checks
explicit settings and mounts, not actual compiled defaults or MPI compatibility.

## Compile once

After making the SIF available, call the following on a node where compilation
is permitted. It can be placed in your own single-task build job:

```bash
bash /absolute/repo/cpp_docker/runner/build.sh /absolute/run/cpppipeline.env
```

The script runs `make clean` and `make` once inside the image. Rebuild after C++
source/header changes. Do not prefix this build script with a many-task `srun`;
all those tasks would otherwise clean and compile the same source tree.

## Run from your own job script

Request the partition, account, task count, memory and time in your own Slurm
script. Request at least two tasks for the distributed pipeline and add this
single invocation after your module setup:

```bash
bash /absolute/repo/cpp_docker/runner/run.sh /absolute/run/cpppipeline.env
```

The wrapper executes **one `srun`**. Do not put another `srun` around it. Slurm
controls the task count and placement; Docker's JSON `ranks` is not used here.
Each task enters the container and executes:

```text
/workspace/src_pipe/Fourier_Quad_Pipe --config /fq-config/pipeline.ini
```

The Bash runner does not inspect job IDs, reject invocation outside an
allocation, submit jobs, infer a queue, run a smoke test or automatically run
`--check`. Slurm and the container runtime report their own launch errors.
Follow your site's rules about where `srun` is permitted.

After changing JSON, rerun `configure` to refresh the env snapshot and INI paths.
Scientific edits outside the managed INI block are retained. Run `--check`
explicitly before production use. Generated env files are trusted Bash input;
do not source another person's unreviewed env file.
