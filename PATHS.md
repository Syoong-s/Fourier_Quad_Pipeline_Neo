# Container path reference

[Home](README.md) · [Docker](DOCKER.md) · [HPC containers](HPC.md)

C++ does not know host-side bind sources. JSON connects two independent pieces:
`mounts` makes host files visible at container destinations; `paths` writes
those destinations into the generated C++ INI. The generator creates the HPC
env from the same mount list. There is no fixed limit on the number of entries.

## Mount once, reuse below it

```json
{
  "mounts": [
    {"host": "/shared/observations", "container": "/observations", "mode": "ro"},
    {"host": "/shared/results", "container": "/results", "mode": "rw", "create": true}
  ],
  "paths": {
    "init.science_root": "/observations/science",
    "init.dq_root": "/observations/dq",
    "init.output_root": "/results"
  }
}
```

This is a fragment to merge into a complete project JSON. Science and DQ share
one parent mount; generated dataset subdirectories share the output mount.
Mount only what your selected workflow needs.

| Mount field | Rule |
| --- | --- |
| `host` | Host directory or ordinary file; relative to the JSON directory when not absolute. |
| `container` | Absolute destination. Paths containing spaces are supported; colons, commas and newlines are not. |
| `mode` | `ro` by default; use `rw` for outputs or mutable data. |
| `create` | Optional Boolean, default false. Only `rw` directories may be created automatically. |

Input paths must exist. Extra mounts need not correspond to a dedicated INI
key; they may serve files named by an exposure list or another mounted file.
Mount an output file's **parent directory** when the program needs to create,
rename or replace it. Read-only input files can be mounted individually.

## Every program path interface

The table covers all runtime file/directory interfaces in both variants. `null`
or an omitted JSON key leaves the INI entry commented, except for the linked
producer/consumer behavior described below. An explicit `""` writes an active
empty string, preserving that field's existing C++ empty-value semantics.

| JSON `paths` key | Use | Notes |
| --- | --- | --- |
| `init.science_root` | Science archive input | Readable directory. |
| `init.dq_root` | DQ archive input | Required when the selected workflow reads DQ. |
| `init.output_root` | Dataset and initialization output | Writable root; includes Science/DQ extractions, stamps, results, lists and astrometry. |
| `astrocat.input_directory` | Raw Gaia input | Flat directory, not recursive. |
| `astrocat.output_directory` | Gaia tile publication | Writable when catalog preparation is enabled. |
| `lensing.astrometry_cat` | Main-phase Gaia tiles | Usually the same mount as Gaia publication. |
| `extcat.input_directory` | Raw source-catalog input | Recursive behavior is a scientific INI setting. |
| `extcat.output_directory` | Source tile publication | Also the main input in Standard; independent in Lite. |
| `lensing.source_cat` | Main source tiles | Alias in Standard; independent input in Lite. |
| `lensing.flat_path` | Super-flat input | Standard only. |
| `lensing.psf_path` | External PSF input | Standard only. |
| `process.expo_list` | Existing top-level exposure list | Input file, not directory; empty derives a per-dataset list. |
| `process.rearr_output_directory` | Rearranged catalog output | Absolute directory or relative name such as `baked`. |
| `process.rearr_output_base_directory` | Base for relative rearrangement output | Empty uses the dataset root derived from input images. |
| `process.rearranged_expo_list_filename` | Published rearranged list name | Filename; an absolute value is a file path, not a directory. |
| `process.rearranged_expo_list_directory` | Published rearranged list directory | Empty uses the input-list parent. |
| `process.fd_expo_list` | FD input list | Input file; empty derives the published rearranged list. |
| `process.fd_output_directory` | FD output | Absolute directory or relative name such as `fdout`. |
| `process.fd_output_base_directory` | Base for relative FD output | Empty uses the dataset root. |

Relative output names do not need their own mount. Do not assign a single fixed
rearrangement/FD destination to multiple datasets unless that sharing is
intentional: per-dataset defaults avoid output collisions.

## Linked paths and aliases

| Relationship | Generator behavior |
| --- | --- |
| Gaia publication ↔ main Gaia input | If only one path is specified, write both INI paths with that value. Two explicit values remain independent. |
| Standard source publication ↔ main source input | One effective C++ field. Conflicting explicit values are rejected; only the authoritative `[extcat].output_directory` key is emitted. |
| Lite source publication ↔ main source input | Independent C++ fields. If only one path is specified, initialize both to it; two explicit values remain independent. |
| Dataset product subdirectories | Derived by C++ under the dataset root; no individual mount required. |

These conveniences belong to the generator; they do not change native C++
configuration semantics. Shared catalog mounts can be read-only for
measurement-only runs and must be writable when their preparation phase runs.

## Additional directories and files

Append a mount to the JSON list:

```json
{"host": "/shared/extra catalog", "container": "/extra/catalog", "mode": "ro"}
```

Or add a read-only list file and refer to it through `process.expo_list`:

```json
{"host": "/shared/lists/expo_gband.list", "container": "/extra/expo.list", "mode": "ro"}
```

Rerun `configure` after changing JSON. Absolute paths inside list files are not
rewritten: mount their referenced dataset paths too. Symbolic-link targets must
also be visible inside the container.

Repeated destinations are rejected. Explicit parent/child mounts are supported;
the deepest destination governs the visible contents and permissions. Be
careful that a child mount hides that part of its parent. System/toolchain,
source and generated-configuration destinations are reserved by the launcher.

## Native runs

None of this JSON is required for a native C++ run. Use host filesystem paths
in your own INI, or configure the headers and rebuild. See the
[C++ guide](CPP_GUIDE.md#configure).
