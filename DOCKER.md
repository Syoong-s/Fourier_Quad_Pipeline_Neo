# Docker

[Home](README.md) · [C++ user guide](CPP_GUIDE.md) · [Container paths](PATHS.md)

Docker is optional. It supplies the C++ toolchain; it does not replace data
preparation or select scientific parameters. The current image targets
`linux/amd64`. Use Linux or a Linux environment such as WSL2 with Python 3.8+
and a working Docker installation. The launcher uses only the Python standard
library and does not require Docker Compose.

## 1. Configure host paths

From the repository root:

```bash
python3 cpp_docker/fqctl.py init
```

Edit `cpp_docker/project.json`. The important fields are:

| Field | Meaning |
| --- | --- |
| `variant` | `standard` or `lite`; must match `source`. |
| `source` | Directory containing the selected C++ Makefile. |
| `image` | Toolchain image reference. Prefer a fixed tag or digest for production. |
| `ranks` | Docker MPI process count; at least 2. |
| `build_jobs` | Concurrent compiler jobs, not MPI ranks. |
| `run_directory` | Generated INI/env directory, separate from the source tree. |
| `mounts` | Host files/directories mapped to container locations. |
| `paths` | Container-side paths assigned to C++ INI keys. |
| `environment` | Explicit container environment; thread counts default to 1. |
| `hpc` | SIF and runtime settings used only for the generated HPC env. |

Relative host paths are resolved against the JSON file's directory, not the
shell's current directory. `project.example.json` shows archive/output mounts;
`project.full.example.json` lists all path interfaces. Neither example is a
complete scientific run configuration. Remove paths and mounts not needed for
the phases you choose. See [Container paths](PATHS.md) for shared mounts and
arbitrary additional inputs/outputs.

## 2. Generate the run configuration

```bash
python3 cpp_docker/fqctl.py configure
```

This generates both files under `run_directory`:

| File | Edit policy |
| --- | --- |
| `pipeline.ini` | Paths come from JSON. Fill/uncomment scientific settings outside the managed-path block. |
| `cpppipeline.env` | Generated HPC launch environment; edit JSON instead of this file. |

Only declared paths are active initially. Other parameters are commented blank
entries such as `# run_process_main =`, so compiled defaults remain effective.
**Standard still defaults to init, main, rearrangement and FD; Lite to init and
main.** Configure the catalogs those phases require, or explicitly disable the
unneeded phases. Review `datasets`, `contains`, geometry, stages and output
policy using the [full reference](CPP_PIPELINE_PARAMETERS.md).

Running `configure` again updates the managed paths, retains your scientific
settings, and saves a backup of the preceding INI. Change paths in JSON, not
inside its managed INI block. A manually created native INI is independent of
this workflow; use the C++ guide to run it directly.

## 3. Check, obtain the image and compile

```bash
python3 cpp_docker/fqctl.py check
python3 cpp_docker/fqctl.py pull
python3 cpp_docker/fqctl.py build
```

`check` inspects JSON, explicit INI settings and declared mounts without starting
Docker or MPI. It does not inspect the effective compiled defaults, all paths
inside exposure lists, FITS content or scientific suitability.

`build` cleans and rebuilds the selected C++ source in the container. Do not
build into one source tree concurrently with another build or a running job.
Header/source changes require another build; scientific INI edits do not.

## Run

```bash
python3 cpp_docker/fqctl.py run
```

The launcher uses `mpirun -np <ranks>` and passes exactly one application option:
`--config /fq-config/pipeline.ini`. This read-only mount refers to the generated
INI, never an INI inside the source directory. The source is read-only during
execution; writable data/output mounts are explicit. Outputs use your host
UID/GID on ordinary Linux bind mounts.

Every command accepts `--json /path/project.json`. For example:

```bash
python3 cpp_docker/fqctl.py run --json /shared/project/project.json
```

## Troubleshooting

| Message or symptom | Action |
| --- | --- |
| JSON changed | Run `configure`, then `check`. Scientific settings are retained. |
| Mount source missing | Correct the host path. Only a `rw` directory with `create: true` is created automatically. |
| Program reads an old site path | Override that path or change the compiled default and rebuild; commented INI entries inherit defaults. |
| Permission denied | Check host permissions, read/write mode and access from the Docker daemon host. |
| Missing executable | Run `build` after selecting the source variant. |
| MPI resource error | Choose a rank count supported by the CPUs assigned to Docker. |
| Image unavailable | Select an accessible compatible toolchain image in JSON; registry access is separate from configuration validation. |
