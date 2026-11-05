"""INI key catalogue; shared by template generation and configuration checks."""
from collections import OrderedDict

# kind: path, relative_path, bool, int, positive, nonnegative, float,
#       text, tokens, datasets, columns, or a | separated set of strings.
COMMON = OrderedDict([
    ('process', OrderedDict([
        ('run_process_astrocat', ('bool', 'Repartition raw Gaia catalogs into one-degree tiles.')),
        ('run_process_extcat', ('bool', 'Repartition raw external source catalogs.')),
        ('run_process_init', ('bool', 'Extract Science/DQ archives and publish datasets and exposure lists.')),
        ('run_process_main', ('bool', 'Run the numerical measurement phase.')),
        ('run_process_rearr', ('bool', 'Spatially rearrange the main-phase catalogs.')),
        ('run_process_fd', ('bool', 'Run scientific field-distortion measurement; not a container diagnostic.')),
        ('expo_list', ('path', 'Existing exposure-list file for downstream-only runs; empty derives from output_root.')),
        ('rearr_output_directory', ('relative_path', 'Rearranged catalog directory; relative to rearr_output_base_directory or the dataset root.')),
        ('rearr_output_base_directory', ('path', 'Rearranged output base; empty uses the dataset root.')),
        ('rearranged_expo_list_filename', ('relative_path', 'Rearranged exposure-list filename, not a directory.')),
        ('rearranged_expo_list_directory', ('path', 'Rearranged list directory; empty uses the input-list parent.')),
        ('fd_expo_list', ('path', 'FD exposure-list file; empty uses the rearranged list.')),
        ('fd_output_directory', ('relative_path', 'FD output directory name or absolute directory.')),
        ('fd_output_base_directory', ('path', 'FD output base; empty uses the dataset root.')),
    ])),
    ('astrocat', OrderedDict([
        ('input_directory', ('path', 'Flat directory of raw Gaia catalogs beginning with RA/Dec columns.')),
        ('output_directory', ('path', 'Published Gaia tiles; usually share one mount with lensing.astrometry_cat.')),
        ('add_header', ('bool', 'true: input starts with data; false: skip exactly one input header line.')),
        ('existing_policy', ('fail|overwrite', 'Existing Gaia tile policy.')),
    ])),
    ('extcat', OrderedDict([
        ('input_directory', ('path', 'Raw external source-catalog directory.')),
        ('output_directory', ('path', 'Published external source-catalog tiles.')),
        ('filename_tokens', ('tokens', 'Comma-separated basename tokens (OR); empty selects all files.')),
        ('recursive', ('bool', 'Search raw catalog subdirectories recursively.')),
        ('delimiter', ('auto|whitespace|comma|tab', 'Input column delimiter.')),
        ('header_mode', ('auto|present|absent', 'Input header handling.')),
        ('malformed_policy', ('fail|skip', 'Malformed-row policy.')),
        ('existing_policy', ('fail|overwrite', 'Existing external tile policy.')),
        ('chunk_mib', ('positive', 'Positive MPI byte-range task size in MiB.')),
        ('total_columns', ('positive', 'Positive external catalog width for pass-through mode.')),
        ('use_explicit_columns', ('bool', 'Set true to use the ordered input_columns projection.')),
        ('input_columns', ('columns', 'Comma-separated one-based input column indices in output order.')),
        ('use_explicit_coordinate_columns', ('bool', 'Set true to use explicit RA/Dec positions rather than header discovery.')),
        ('ra_column', ('positive', 'One-based RA column in the raw input.')),
        ('dec_column', ('positive', 'One-based Dec column in the raw input.')),
        ('mag_g_column', ('nonnegative', 'One-based g-magnitude column; 0 marks the band absent.')),
        ('mag_r_column', ('nonnegative', 'One-based r-magnitude column; 0 marks the band absent.')),
        ('mag_i_column', ('nonnegative', 'One-based i-magnitude column; 0 marks the band absent.')),
        ('mag_z_column', ('nonnegative', 'One-based z-magnitude column; 0 marks the band absent.')),
        ('mag_y_column', ('nonnegative', 'One-based y-magnitude column; 0 marks the band absent.')),
        ('zp_column', ('positive', 'One-based zp/photo-z column in the raw input; not a magnitude zero point.')),
    ])),
    ('init', OrderedDict([
        ('science_root', ('path', 'Raw Science FITS/FZ archive root.')),
        ('dq_root', ('path', 'Raw DQ FITS/FZ archive root.')),
        ('output_root', ('path', 'Dataset output root; science, dqmask, stamps, result, astrometry and expolists share this mount.')),
        ('datasets', ('datasets', 'Comma-separated TARGET:PREFIX pairs, e.g. gband:c4d_; replaces the old dataset/target/prefix options.')),
        ('contains', ('tokens', 'Comma-separated archive basename tokens (OR); an active empty value disables filtering.')),
        ('existing', ('fail|resume|overwrite', 'Existing dataset policy.')),
        ('f77_max_path', ('nonnegative', 'Maximum generated path length; 0 disables the guard.')),
    ])),
])
LENSING_COMMON = OrderedDict([
    ('process_stage', ('positive', 'Product of enabled stage primes. Factor 23 requires factor 19.')),
    ('astrometry_cat', ('path', 'Gaia tile directory read by the main phase.')),
    ('source_cat', ('path', 'Source tiles: aliases extcat.output_directory in Standard; independent in Lite.')),
    ('ccd_split', ('int', 'CCD partition: 1 for the whole CCD, 2 for amplifier halves.')),
    ('gal_smooth', ('nonnegative', 'Galaxy-stamp smoothing selector; non-negative.')),
    ('star_smooth', ('nonnegative', 'Star-stamp smoothing selector; non-negative.')),
    ('pixel_size', ('float', 'Pixel scale in arcsec/pixel; positive.')),
    ('nmax_chip', ('positive', 'Positive CCD count.')),
    ('chipnx', ('positive', 'Physical CCD width in pixels; positive.')),
    ('chipny', ('positive', 'Physical CCD height in pixels; positive.')),
])
STANDARD_ONLY = OrderedDict([
    ('astrometry_trivial', ('int', 'Astrometry: 0 for Gaia, 1 for identity mapping (Standard only).')),
    ('astrometry_cat_type', ('int', 'Gaia layout: 1 for legacy large tiles, 2 for one-degree tiles (Standard runtime key).')),
    ('include_flat', ('int', 'Super-flat correction: 0 off, 1 on (Standard only).')),
    ('include_mask', ('int', 'Mask mode: 0, 1, 2 or 3 (Standard only).')),
    ('flat_path', ('path', 'Directory containing super-flat FITS images (Standard only).')),
    ('psf_path', ('path', 'External PSF FITS root (Standard only).')),
    ('ext_cat', ('int', 'External sources: 0 off, 1 on; FD requires 1 (Standard only).')),
    ('ext_psf', ('int', 'PSF input: 0 for exposure stars, 1 for external PSF (Standard only).')),
    ('deblending', ('int', 'Deblending: 0 off, 1 on (Standard only).')),
    ('psf_type', ('int', 'PSF model: 1 local polynomial, 2 hybrid (Standard only).')),
    ('psf_ms', ('int', 'PCA/multi-scale PSF reconstruction: 0 off, 1 on (Standard only).')),
])

def sections(variant):
    if variant not in ('standard', 'lite'):
        raise ValueError('variant must be standard or lite')
    result = OrderedDict((s, OrderedDict(values)) for s, values in COMMON.items())
    result['lensing'] = OrderedDict(LENSING_COMMON)
    if variant == 'standard':
        result['lensing'].update(STANDARD_ONLY)
    return result


def path_keys(variant):
    return {s + '.' + k: kind for s, values in sections(variant).items()
            for k, (kind, _) in values.items() if kind in ('path', 'relative_path')}


# Always writable when explicitly mapped. Catalog outputs can instead be readonly
# when their preparation stage is not selected; they are also main-stage inputs.
OUTPUT_PATHS = {
    'init.output_root', 'process.rearr_output_directory',
    'process.rearr_output_base_directory', 'process.rearranged_expo_list_directory',
    'process.fd_output_directory', 'process.fd_output_base_directory',
}
FILE_PATHS = {'process.expo_list', 'process.fd_expo_list'}

# An absolute rearranged-list filename is a file path, not a directory.
FILE_PATHS.add('process.rearranged_expo_list_filename')
OUTPUT_PATHS.add('process.rearranged_expo_list_filename')
