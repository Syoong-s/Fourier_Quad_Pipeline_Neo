# Maintenance notes

## Scope

This change is based on `9925f4f4f19b16344274f8c4f85d3fb785598423`. It updates
configuration interfaces, container launch helpers and documentation. It does
not replace `src/process_*`, scientific configuration headers, Makefiles or the
FITS decompression correction in that commit.

| File | Responsibility |
| --- | --- |
| `cpp_docker/schema.py` | INI key catalogue, type hints and path classification. |
| `cpp_docker/fqctl.py` | JSON validation, path linking, INI/env generation and Docker argv. |
| `cpp_docker/runner/common.sh` | Trusted env loading and explicitly requested checks. |
| `cpp_docker/runner/run.sh` | One `srun`; no allocation/step guard or automatic config check. |
| `cpp_docker/runner/runtime.sh` | One container per rank, its actual PMI environment and explicit binds. |
| `cpp_docker/runner/build.sh` | One containerized build without scheduler-specific gating. |
| `cpp_*/src/RuntimeConfig.cpp` | Compiled defaults, transactional INI parsing, config-only CLI. |
| `cpp_*/main.cpp` | Existing MPI configuration broadcast and scientific dispatch. |

The internal legacy `help_requested` output remains false for compatibility with
existing internal signatures; no help CLI is accepted. The former
`parseRuntimeCommandLine` entry only validates the config-only grammar. Explicit
INI `process.expo_list` updates the previous explicit-list state.

## Generation contract

INI paths occupy a marked block; science settings outside it are preserved on
regeneration. Unknown or duplicate active keys are rejected by the launcher.
Repeated section headings are supported by C++ and permit separate path/science
blocks. Commented blanks are intentionally distinct from active empty values.

Generation uses an advisory lock and atomic replacement for each file, retaining
an INI backup. An identifier ties the INI/env pair together; the JSON digest
helps the explicit configuration check detect stale generated settings. These
are runtime project consistency controls, not an upgrade installer.

Mount coverage uses the deepest matching destination. JSON can share a parent
mount or declare explicit children. System, source and configuration destinations
are reserved. Standard has one source-catalog runtime field; Lite and the Gaia
producer/consumer pairs remain independently configurable in C++.

A configuration check covers only explicit settings and declared mounts. It
cannot infer all compiled constants, validate every file referenced in lists,
validate scientific results or certify MPI transport compatibility.

## Runtime boundary

Docker argv is built as an argument list, not an evaluated shell string. Bind
sources must exist unless explicitly declared creatable output directories.
Runs use current UID/GID and supplemental groups, source/config read-only, and
explicit writable data mounts.

Each Slurm task enters `runtime.sh` before its own `SLURM_*`, `PMI_*` and
`PMI2_*` values are forwarded. The script retains environment-prefix cleanup,
`--cleanenv`, `--no-eval` and explicit bindings for predictable container MPI
startup. This is not a sandbox for malicious env files. The runner sources
trusted generated Bash and leaves resource decisions to the user's scheduler
script. No `SLURM_JOB_ID`, nested-step or task-count guard is applied.

## Regression checks

The tests are maintenance checks, not smoke-test or dry-run commands exposed to
pipeline users:

```bash
python3 -m unittest discover -s cpp_docker/tests -v
```

In a full repository checkout, compile the configuration tests with real headers:

```bash
set -euo pipefail
for variant in Standard Lite; do
  g++ -std=c++17 -Wall -Wextra -Werror \
    -I"cpp_${variant}" -I"cpp_${variant}/include" -I"cpp_${variant}/config" \
    "cpp_${variant}/src/RuntimeConfig.cpp" \
    "cpp_${variant}/tests/RuntimeConfigTest.cpp" \
    -o "/tmp/fq-config-${variant}"
  "/tmp/fq-config-${variant}" "cpp_${variant}/pipeline.example.ini"
done
```

The CI configuration uses this real-header command. Offline launch tests use
substitute Docker/srun/Apptainer executables and do not demonstrate actual image
startup or multi-node MPI. A separate complete build and representative
scientific dataset are required for numerical acceptance.

## Image provenance

Dockerfile and Pixi lock files are unchanged. The build installs a locked
G++ 12.3.0 scientific stack, OpenMPI 4.1.8 and the Slurm 25.11.2 PMI2 client.
It does not embed pipeline source, observation data, a Slurm controller,
scheduler command set, MUNGE or site configuration.

| Base | Pinned reference |
| --- | --- |
| Pixi | `ghcr.io/prefix-dev/pixi:0.71.1@sha256:cc1401e1553164ad032908e1f6a1dbbaa8c8e067ce4347ae831c43a066bbc66c` |
| Rocky Linux | `docker.io/rockylinux/rockylinux:8.10@sha256:f5529992e67440c1a4ae7788244d4381c6909159a88eacd95b7523ae47ced82e` |

| Source archive | SHA-256 recorded in Dockerfile |
| --- | --- |
| `https://download.schedmd.com/slurm/slurm-25.11.2.tar.bz2` | `9bfd844f746c1268d6e0cb14bda586b548b1c4b12f8c72181397d0cec9ce5d39` |
| `https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.8.tar.bz2` | `466f68e3132a1dc02710cc2011fafced8336d98359fa2dae4dddcfd5719f12a9` |

`.dockerignore` excludes project files/data from the toolchain build context.
The release workflow packages the complete repository so root navigation and
scientific source ship together. Third-party notices remain in the root notice
file; dependencies retain their upstream licenses.

## Local validation of this replacement package

| Check | Result and scope |
| --- | --- |
| Launcher and documentation unit tests | 58 passed; 1 skipped because the replacement-only tree does not contain unchanged scientific headers. |
| Standard configuration parser | 136 assertions passed with isolated fixture headers; not a complete scientific build. |
| Lite configuration parser | 136 assertions passed with isolated fixture headers; not a complete scientific build. |
| Parameter-reference runtime-key coverage | Every Standard/Lite INI key in the launcher catalogue is documented. |
| Actual-header inventory check | Included in CI; requires the full base checkout and is the skipped local test above. |
| Bash, Python, JSON and workflow syntax | Checked locally. |
| Real Docker/SIF/Slurm and scientific results | Not executed; those runtimes and observation data are unavailable in this environment. |

The private parser fixtures and compiled test executables are not distributed
in this replacement package and do not replace configuration headers. Root
Markdown, source/template comments in supplied files, examples and generated
messages were checked for untranslated Chinese text. ZIP extraction was checked
against the supplied file contents. A complete clean-checkout integration and
numerical validation still require the real repository and scientific stack.
