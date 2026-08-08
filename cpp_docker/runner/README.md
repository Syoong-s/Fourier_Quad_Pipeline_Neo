# cpp_Standard on Slurm/Apptainer clusters

This runner consumes one precompiled SIF on any x86_64 Slurm cluster that
advertises the PMI2 launcher plugin. The SIF supplies G++ 12.3.0,
OpenMPI 4.1.8, a PMI2 client, and the scientific stack. Host OpenMPI and host
compiler ABI compatibility are not part of the launch contract.

The process boundary is:

`srun --mpi=pmi2` → `run-apptainer.sh` → `apptainer exec --cleanenv` →
container-linked `Fourier_Quad_Pipe`.

## Site prerequisites

- x86_64 Linux compute nodes;
- Slurm with `pmi2` listed by `srun --mpi=list`;
- Apptainer or Singularity available on compute nodes;
- one shared filesystem visible at identical paths on all allocated nodes;
- routable TCP between allocated nodes for the portable baseline.

The image does not include vendor UCX, OFI, or RDMA providers. Those are
performance extensions that require separate qualification.

## Directory and configuration

A typical shared layout is:

```text
/shared/project/cpppipeline/
├── code/
├── runner/
├── images/
├── apptainer-cache/
├── apptainer-tmp/
├── scratch/
└── data/
    ├── Science/
    ├── DQMask/
    ├── AstroDir/
    ├── ExtSrcDir/
    ├── FlatDir/
    └── DataProcess/
```

Copy this complete runner directory and create the trusted configuration:

```text
cp cpppipeline.env.example cpppipeline.env
```

Edit every host path. `CPP_SIF`, source, data, scratch, and the runner must be
visible from every allocated node. Container catalogue destinations must match
the strings compiled into
`/workspace/src_pipe/include/process_main/LensingConfig.hpp`.

`SCIENCE_ROOT_HOST` and `DQ_ROOT_HOST` are optional read-only archive binds.
Set them only when running process_init or process_main. Their container
destinations are `SCIENCE_ROOT_CONTAINER` and `DQ_ROOT_CONTAINER`; use exactly
those values for the executable's `--science-root` and `--dq-root` arguments
when active. Four additional optional mounts (`EXTCAT_INPUT`, `REARR_OUTPUT`,
`EXPOLIST_DIR`, `FD_OUTPUT`) follow the same conditional pattern: set `*_HOST`
to bind, leave empty to skip. Generated files remain under the writable
`PROCESS_DATA_CONTAINER`.

`HPC_MODULES`, `HPC_EXTRA_BINDS`, `HPC_PASSTHROUGH_ENV`,
`HPC_CONTAINER_ENV`, and `SRUN_ARGS` are Bash indexed arrays. Modules may make
Apptainer or Slurm available, but must not inject a host MPI into the
application environment. The runner uses `--cleanenv`, forwards every
`SLURM_*`, `PMI_*`, and `PMI2_*` value created by Slurm, and forwards only
the explicitly configured extra environment.

## Acquire the SIF

For a reviewed Docker archive, set `CPP_DOCKER_ARCHIVE`, `CPP_SIF`,
`APPTAINER_CACHE_DIR`, and `APPTAINER_TMP_DIR`, then submit:

```text
sbatch build-sif.slurm
```

For a registry image, set a digest-pinned `OCI_IMAGE_URI` and run:

```text
bash pull-sif.sh
```

Both paths build or pull to a temporary file, atomically rename the finished
SIF, create `${CPP_SIF}.sha256`, and refuse existing outputs. Copy the first
field of that sidecar to `CPP_SIF_SHA256_EXPECTED` before production use.

## Validation order

Run the singleton image and bind check:

```text
bash run-apptainer.sh --check
```

Compile the full pipeline once on a compute node:

```text
sbatch compile-pipeline.slurm
```

Run a one-node, two-rank smoke test:

```text
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
```

Then run the default multi-node smoke test:

```text
sbatch mpi-smoke-test.slurm
```

The smoke job compiles self-contained MPI and scientific-stack probes in a
unique directory below `HPC_SHARED_SCRATCH_HOST`. It requires at least two
ranks and removes only its own temporary directory.

After reviewing catalogue paths and the exposure list, launch the pipeline:

```text
sbatch cpppipeline.slurm
```

Without script arguments, the executable receives
`${PROCESS_DATA_CONTAINER}/expo_list.list`. Additional arguments after the
script name are passed to `Fourier_Quad_Pipe`.

For a chained initialization and numerical run, pass the named workflow options,
for example `--run-init true --run-main true --science-root
/data/archive/science --dq-root /data/archive/dqmask --output-root
/data/DataProcess --dataset g2013:c4d_13 --dataset g2014:c4d_14 --contains v1`.
Repeated datasets run sequentially; repeated contains tokens use OR matching.

## Scheduler templates

The `#SBATCH` resources are conservative templates. Override node, task,
CPU, memory, partition, account, and time values together according to local
policy. If a site requires centralized logs, pass absolute `--output` and
`--error` paths to every `sbatch` command.

Do not compile the same source copy concurrently. Source and processing binds
are writable; catalogue and calibration binds are read-only. Science, DQMask,
and four pipeline output mounts are optional (bound only when `*_HOST` is set).
