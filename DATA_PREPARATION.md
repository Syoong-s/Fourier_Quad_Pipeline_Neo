# Data preparation

[Home](README.md) · [C++ user guide](CPP_GUIDE.md) · [Parameter reference](CPP_PIPELINE_PARAMETERS.md)

Prepare these inputs before running the scientific pipeline. Paths in a native
run are ordinary filesystem paths. Paths used inside a container must be visible
at the corresponding container locations.

## Required inputs

| Input | Purpose | Minimum requirement |
| --- | --- | --- |
| Science images | Exposures used for source detection, shape measurement, and weak-lensing processing. | Supported Science FITS/FZ data matching the configured archive, exposure, and CCD naming rules. |
| Gaia catalog | Sky-coordinate reference for source matching and astrometric calibration. | Covers the Science footprint and follows the configured Type 1 or Type 2 tile layout. |
| External source catalog | Supplies sky positions, `zp`, and photometry for matching and downstream processing. | Canonical one-degree tiles containing `ra`, `dec`, `zp`, and at least one observed-band magnitude. |
| DQ masks (configuration-dependent) | Mark bad, saturated, defective, or otherwise invalid pixels. | Required by Lite and by every Standard configuration that reads DQ data. |

These are the normal measurement inputs. A catalog-preparation-only run needs
only its raw catalog input and output directory. Standard branches that bypass
an input must also disable every active stage that otherwise consumes it.

## Science images

Science images are the primary scientific exposures, not calibration catalogs or
an output directory. They must use a FITS/FZ format supported by the selected
variant, be readable through its FITS logic, and follow the exposure/CCD naming
and list conventions below.

The initializer searches the configured `[init].science_root` and `dq_root`.
`datasets` contains `TARGET:PREFIX` pairs; `TARGET` names the output dataset and
`PREFIX` selects archive basenames. `contains` supplies additional case-sensitive
basename tokens using OR matching. An active `contains =` disables that token
filter; a commented entry keeps the compiled filter, which is `v1` by default.

| Convention | Default | Where to configure |
| --- | --- | --- |
| Archive suffix | `.fits.fz` | `InitConfig::ARCHIVE_SUFFIX` |
| CCD identity keyword | `CCDNUM` | `InitConfig::CCDNUM_KEYWORD` |
| DQ stem replacement | Every `ood` becomes `ooi` | `DQ_STEM_REPLACE_FROM` / `DQ_STEM_REPLACE_TO` |
| Dataset | `gband:c4d_` | `[init].datasets` |
| Basename token | `v1` | `[init].contains` |

The archive/FITS naming constants are in Standard's `config/InitConfig.hpp` and
Lite's `Initialize.hpp`; changing them requires rebuilding. Extracted chip
products are identified by physical `CCDNUM`, not by their HDU position. Use
matching Science/DQ exposures and detector identities. Do not mix products from
older HDU-order naming with the current output layout.

## Gaia catalog

The Gaia catalog supplies accurate RA/Dec reference positions for object
matching and astrometric calibration. It must cover the actual Science-image
footprint and be stored directly under the configured `ASTROMETRY_CAT`
directory. Every consumed tile starts with one header line; subsequent rows use
the first two numeric fields as RA and Dec. Rows may be comma- or
whitespace-separated and additional fields are ignored.

Standard selects its filename layout with `LensingConfig::AstroCatType` or
`[lensing].astrometry_cat_type`. Type `1` is the legacy large-tile layout below.
Type `2` uses one-degree tiles and may be generated from raw flat-directory
files by `process_astrocat`. **Lite always uses the Type-2 one-degree layout**;
it has no `AstroCatType` symbol or `astrometry_cat_type` INI key.

### Type 1: legacy large tiles

- `|Dec| < 80°`: `gaia_<p|m><D>_<RR>.cat`, where
  `D = floor(|Dec| / 10) + 1` (1–8) and `RR = floor(RA / 10)` (00–35,
  zero-padded).
- `|Dec| >= 80°`: `gaia_<p|m>9.cat`, without an RA suffix.
- `p` denotes nonnegative Dec; `m` denotes negative Dec. Each file contains one
  header line.

| Filename | Coverage |
| --- | --- |
| `gaia_p1_00.cat` | `0° <= RA < 10°`, `0° <= Dec < 10°` |
| `gaia_m3_12.cat` | `120° <= RA < 130°`, `-30° < Dec <= -20°` |
| `gaia_p9.cat` | `0° <= RA < 360°`, `80° <= Dec <= 90°` |

For Type 1, to ensure that stars can still be selected for exposures located at
the edges of the 10° × 10° grid, expand the upper and lower Dec coverage limits
of a single star catalog by 2°. Expand the upper and lower RA limits by 2°, 4°,
and 6° within the 0°, 30°, and 60° ranges, respectively, as in the existing
catalog-preparation convention.

### Type 2: one-degree tiles

Type 2 files use
`<ASTROMETRY_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`. The default
prefix in `config/pathconfig.hpp` is `astra_`; it does not include `RA_`.
RA bounds use three digits and Dec bounds use signed `p`/`m` two-digit values.
For example, `astra_RA_123_124_Dec_m05_m04.dat` covers
`123° <= RA < 124°`, `-5° <= Dec < -4°`.

`process_astrocat` reads only direct regular files from
`[astrocat].input_directory`, produces these Type-2 tiles, and removes
exact/one-ULP coordinate duplicates. Its `output_directory` controls only
publication. To consume the result, set `[lensing].astrometry_cat` to that
directory and use the appropriate layout:

| Variant | Select Type 2 |
| --- | --- |
| Standard | Set `[lensing].astrometry_cat_type = 2` in the INI. |
| Lite | Already fixed to Type 2; set the consumer path only. There is no layout selector to edit. |

Set `[astrocat].add_header = true` when each raw file starts with data, or
`false` to skip exactly one input header line. Output tiles always have a
header. Use finite coordinates with `0 <= RA <= 360` and `-90 <= Dec <= 90`;
RA 360 is stored as 0 and Dec 90 belongs to the final Dec tile.

## External source catalog

The minimum schema is intentionally survey- and band-independent:

| Field | Meaning |
| --- | --- |
| `ra` | Right Ascension. |
| `dec` | Declination. |
| `zp` | The catalog `zp` quantity consumed by the selected Pipeline configuration; the FD cuts treat this as photometric redshift, not a magnitude zero point. |
| One observed-band magnitude | A magnitude in any one band used by the selected analysis. |

Additional colors, redshifts, object classes, shapes, and flags may be retained,
but they are not part of this minimum input contract. Configure actual column
positions, delimiter, header handling, and projection in `[extcat]` or in
`config/ExtCatConfig.hpp` before rebuilding.

### Tile filenames

| Component | Convention |
| --- | --- |
| Tile size | 1° × 1° |
| Filename | `<SOURCE_CAT_TILE_PREFIX>RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat` |
| Default prefix | `extern_` in `config/pathconfig.hpp`, excluding the fixed `RA_` token |
| RA boundaries | Three digits |
| Dec boundaries | `p` or `m` followed by a two-digit absolute value |
| Upper boundary | One degree above the lower boundary |
| Header | Exactly one header line before data |

For example, `extern_RA_123_124_Dec_m05_m04.dat` covers
`123° <= RA < 124°` and `-5° <= Dec < -4°`.

Use `process_extcat` to convert raw downloads into these canonical pipeline
tiles. Raw downloader filenames are source artifacts and are not a second C++
naming convention. The output must not equal or be nested below the raw input.

### Column mapping

Column numbers in `[extcat]` are **one-based indices of the raw input**. With
explicit projection, `input_columns` specifies the output order; retain RA, Dec,
`zp` and every configured, nonzero magnitude field. Set
`use_explicit_columns = true` to enable that projection and
`use_explicit_coordinate_columns = true` when supplying explicit RA/Dec
positions instead of header discovery.

A magnitude column of `0` marks that band absent. FD requires at least one
available band and selects the first available in **i → z → r → g → y** order.
Do not hand-edit downstream column offsets to match a new raw catalog.

In Standard, `[extcat].output_directory` is also the main pipeline's source-tile
input; `[lensing].source_cat` is an alias for the same value. In Lite, the producer
output and `[lensing].source_cat` are separate fields: set both when consuming
newly generated tiles.

## DQ masks

DQ masks identify pixels that must not participate in scientific measurements,
including bad pixels, saturation, detector defects, and other invalid regions.
They are optional only when the chosen configuration disables DQ access. Standard
users can select the relevant mask mode with `include_Mask`; current C++ Lite is
fixed to per-chip DQ masks, so Lite runs must provide them. If DQ masks are
omitted, verify that no configured path or active branch still reads them.

Use masks compatible with the Science images' geometry, exposure names and
`CCDNUM` identities. Check the selected initialization path as well as the main
measurement mask mode before deciding that a DQ archive is unnecessary.

## Existing exposure lists

For downstream-only processing, set `[process].expo_list` to an existing top-level
list and disable `run_process_init`. Each nonblank record contains a per-exposure
chip-list path; a trailing legacy chip count is accepted. Each chip list contains
the Science image paths. The current reader supports quoted paths.

The main/rearrangement phases derive the dataset root from image paths, not
solely from `[init].output_root`. Keep the established dataset layout. Changing
an output-root setting alone does not relocate an existing dataset or rewrite
its lists. Inside containers, list entries must name mounted container paths.

For multiple downstream-only datasets, leave `expo_list` empty and let the
program derive `output_root/expo_<target>.list` for each dataset; one explicit
list cannot serve multiple main/rearrangement datasets.

## Before the first measurement run

Check archive naming and detector identities, catalog sky coverage, tile names,
column mapping, enabled phases and output locations. Use a small representative
set of exposures before starting a full survey run. Input directories should be
readable; every selected output location should be writable.
