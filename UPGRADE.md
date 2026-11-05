# Upgrade from `9925f4f`

[Home](README.md) · [C++ user guide](CPP_GUIDE.md)

Base commit: `9925f4f4f19b16344274f8c4f85d3fb785598423`.
This is a replacement-file package for that checkout, not a complete source
archive and not a patch against the previously supplied ZIP.

## Replace files

1. Back up local source/header changes, INI files, runner configuration and data.
2. Copy the package contents into the repository root, **merging directories**
   and replacing same-name files. Include the hidden `.github` directory.
3. Delete the obsolete files listed in `REMOVE.txt` when they exist. Plain ZIP
   extraction does not delete them.
4. Rebuild the chosen C++ variant before running it.

**Do not delete `cpp_Standard`, `cpp_Lite` or `cpp_docker` before copying.**
Their unchanged algorithms, configuration headers, Makefiles, Dockerfile and
lock files are required and are deliberately absent from this package. The
FITS decompression changes already in `9925f4f` are not replaced.

No installer, patch-application helper or package-verification script is
required. `REMOVE.txt` is a list to review, not an executable deletion script.

## What remains unchanged

Scientific algorithms and compiled defaults are retained. A no-argument C++
invocation uses only compiled settings; it does not look for an INI. Explicit
empty INI values keep their existing behavior: they are not a general request
to inherit defaults. `process_fd` remains a scientific pipeline phase.

## Move old application options into the INI

The executable now accepts only `--config PATH` or `--config=PATH`. Old help,
scientific-option and positional-list arguments are rejected. MPI launcher
arguments such as `mpirun -np 4` are unaffected.

| Former option | INI replacement |
| --- | --- |
| `--run-astrocat`, `--run-extcat`, `--run-init`, `--run-main`, `--run-rearr`, `--run-fd` | Corresponding `[process].run_process_*` Boolean. |
| `--astrocat-input`, `--astrocat-output`, `--astrocat-add-header`, `--astrocat-existing` | `[astrocat].input_directory`, `output_directory`, `add_header`, `existing_policy`. |
| `--extcat-input`, `--extcat-output`, `--extcat-contains`, `--extcat-recursive` | `[extcat].input_directory`, `output_directory`, `filename_tokens`, `recursive`. |
| `--extcat-delimiter`, `--extcat-header`, `--extcat-malformed`, `--extcat-existing`, `--extcat-chunk-mib` | `[extcat].delimiter`, `header_mode`, `malformed_policy`, `existing_policy`, `chunk_mib`. |
| `--extcat-columns` | `[extcat].input_columns` **and** `use_explicit_columns = true`. |
| `--extcat-ra-column`, `--extcat-dec-column` | `[extcat].ra_column`, `dec_column` **and** `use_explicit_coordinate_columns = true`. |
| `--extcat-zp-column` | `[extcat].zp_column`. |
| `--science-root`, `--dq-root`, `--output-root` | `[init].science_root`, `dq_root`, `output_root`. |
| Repeated `--dataset`, or legacy `--target` and `--prefix` | `[init].datasets = TARGET:PREFIX,...`. |
| Repeated `--contains`, `--existing`, `--f77-max-path` | `[init].contains`, `existing`, `f77_max_path`. |
| `--expo-list` or a bare list path | `[process].expo_list`. |
| `--rearr-output-dir`, `--rearr-output-base`, `--rearr-list-name`, `--rearr-list-dir` | `[process].rearr_output_directory`, `rearr_output_base_directory`, `rearranged_expo_list_filename`, `rearranged_expo_list_directory`. |
| `--fd-expo-list`, `--fd-output-dir`, `--fd-output-base` | `[process].fd_expo_list`, `fd_output_directory`, `fd_output_base_directory`. |

Native users may keep their existing INI and run it directly. The new repository
INI templates have commented blank entries so that they do not silently select
a different scientific configuration.

## Container users

Create a JSON project with the [Docker guide](DOCKER.md), generate its separate
INI/env and review the enabled phases. Generation writes only configured paths;
other fields begin as `# key =`. The generated INI is passed explicitly with
`--config`, separate from any native INI beside the executable.

For HPC, call the new `runner/run.sh` once from your own job script. It delegates
task counts and allocation behavior to `srun`; it has no job-ID/step checks or
automatic configuration-test invocation. Use `run.sh --check` explicitly when
needed. Old Compose files, complete Slurm templates and runtime smoke-test
utilities are removed; their removal does not affect scientific FD analysis.

All user guides are English files at the repository root. Data preparation and
the full C++ parameter reference remain available through `README.md`.
