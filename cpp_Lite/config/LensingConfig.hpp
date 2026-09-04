#ifndef LENSING_CONFIG_HPP
#define LENSING_CONFIG_HPP

#include "Initialize.hpp"
#include "pathconfig.hpp"

#include <cmath>

// ==========================================
// cpp_lite: the following build-time branch selectors of the full pipeline are FROZEN and
// their unused branches have been deleted from the sources.  They no longer exist as
// constants; the code implements exactly one behaviour each:
//   ASTROMETRY_trivial = 0  -> Gaia-based astrometry only
//   include_FLAT       = 0  -> no super-flat multiplication
//   include_Mask       = 2  -> per-chip DQ mask from dirOutput/dqmask
//   ext_cat            = 1  -> external (SOURCE_CAT) source catalogue
//   ext_PSF            = 0  -> PSF measured from the stars of the frame
//   deblending         = 1  -> de-blending always applied
//   PSF_type           = 1  -> local polynomial PSF fit
//   PSF_Ms             = 0  -> no multi-scale / PCA PSF reconstruction
// Runtime-selectable Lite fields also include the Gaia/source paths, pixel
// scale, chip count, and physical chip geometry; see pipeline.example.ini.
// ==========================================
namespace LensingConfig {
    // Compiled camera fallbacks copied into RuntimeConfig before file/CLI overrides.
    inline constexpr int DEFAULT_CHIP_COUNT = Initialize::N_CCD;
    inline constexpr double pixel_size = Initialize::pixel_size;
    inline constexpr double saturation_thresh = Initialize::saturation_thresh;
    inline constexpr int chipnx = Initialize::chipnx;
    inline constexpr int chipny = Initialize::chipny;

    // Stage control parameters
    inline constexpr int PROCESS_stage = Initialize::PROCESS_stage;

    // Split parameters
    inline constexpr int CCD_split = Initialize::CCD_split;
    constexpr int nct = 12;  // Number of background rectangles.
    constexpr int ncx = 3;  // Number of background rectangles along x.

    // PSF selection and configuration
    constexpr int npo = 64;  // Exposure PSF sample count.
    constexpr int nstar_min = npo * 3 / 2;  // Minimum stars for exposure PSF fitting.
    constexpr int npl = 10;  // Local PSF polynomial coefficient count minus one.
    constexpr int nstar_min_local = 16;  // Minimum retained stars for a local fit.

    // ==========================================
    // Configuration: Stage-5 stellar-locus, grouping, and PRESS selection
    // Method: Keep the shared scientific cuts compile-time selectable while
    //         Lite's sole Gaia/local branch remains structurally fixed.
    // ==========================================
    constexpr int psf_exposure_min_candidates = 60;  // Minimum exposure-wide PSF candidates.
    constexpr double psf_count_pilot_clip_sigma = 3.0;  // Robust star-area pilot clipping multiplier.
    constexpr int psf_count_pilot_clip_iterations = 3;  // Robust star-area pilot clipping passes.
    constexpr double psf_count_zero_mad_quantile = 0.05;  // Symmetric lower quantile; upper is one minus this value.
    constexpr double psf_count_hist_range_sigma = 5.0;  // Local star-area histogram half-range in pilot widths.
    constexpr double psf_count_locus_sigma = 4.0;  // Exposure star-area locus sigma window.
    constexpr int psf_count_locus_min_samples = 30;  // Minimum star-area locus samples.
    constexpr double psf_minchi_reference_fraction = 1.0 / 3.0;  // Exposure top-size reference fraction.
    constexpr int psf_minchi_reference_max_per_chip = 5;  // Reference-star cap per chip.
    constexpr double psf_minchi_sigma_cut = 4.0;  // Minimum-chi rejection sigma.
    constexpr double psf_gaia_match_radius_pix = 2.0;  // Gaia match radius in pixels.
    constexpr int psf_gaia_locus_min_matches = 5;  // Minimum Gaia matches for locus support.
    constexpr double psf_pair_chi_valid_peak_fraction = 0.3678794411714423216;  // exp(-1) pair-chi peak threshold.
    constexpr double psf_bad_fraction_valid_peak_fraction = 0.10;  // Bad-pair-fraction peak threshold.
    constexpr double psf_type3_elbow_search_height_fraction = 0.10;  // Elbow candidates must lie below this smoothed main-peak fraction.
    constexpr bool psf_press_rejection_enabled = true;  // Enable optional post-fit PRESS cleanup.
    constexpr double psf_press_sigma_cut = 4.0;  // Standardized PRESS rejection sigma.
    constexpr int psf_press_max_removals = 5;  // Maximum PRESS removals permitted per chip.
    constexpr double psf_loo_min_denom = 1.0e-6;  // Minimum leave-one-out denominator.
    // ==========================================

    // Stamp dimensions
    constexpr int ns = 64;  // Science stamp and Fourier-grid side length.
    constexpr int nsns = ns * ns;  // Pixels in one science stamp.
    constexpr int chip_margin = 8;  // Extra chip-edge extraction margin.
    constexpr int ns_2 = ns / 2;  // Half science-stamp side length.
    constexpr int nl_2 = ns_2 + chip_margin;  // Half expanded extraction side.
    constexpr int nl = nl_2 * 2;  // Full expanded extraction side.
    constexpr int flag_thresh = 3;  // Maximum accepted source extraction flag.
    constexpr int chip_edge_margin = chip_margin;  // Alias used by chip-edge checks.

    constexpr double dz_thresh = 0.1;  // Redshift tolerance for catalog matching.

    // Initial catalog reservation capacities. These values are not hard
    // source or star limits; vectors may grow beyond them as needed.
    // Maximum number of flux-ranked image detections passed to astrometric pattern matching.
    // This is a scientific selection limit, not a catalog-storage capacity limit.
    constexpr int n_user_max = 200;  // Bright detections used for astrometric matching.
    constexpr int ngal_max = 2000;  // Initial galaxy-vector reservation hint.
    constexpr int nstar_max = 1000;  // Initial star-vector reservation hint.
    constexpr int src_npara = 12;  // Source-catalog metadata field count.

    constexpr int npd = 33;  // PU astrometric distortion coefficient count.
    // Target side length for the balanced background blocks.
    constexpr int blocksize = 200;  // Target background block side length.
    constexpr int bg_rough_grid_x = 32;  // Rough background grid columns.
    constexpr int bg_rough_grid_y = 32;  // Rough background grid rows.
    constexpr int bg_min_block_pixels = 1000;  // Minimum pixels in a background block.
    constexpr int bg_min_clipped_pixels = 200;  // Minimum pixels after block clipping.
    constexpr double bg_min_valid_frac = 0.25;  // Minimum valid fraction per background block.
    constexpr double bg_clip_low = 4.0;  // Lower background clipping sigma.
    constexpr double bg_clip_high = 2.5;  // Upper background clipping sigma.
    constexpr double bg_fit_clip_sigma = 3.0;  // Background-plane fit clipping sigma.
    constexpr int bg_fit_max_iter = 4;  // Maximum background-plane clipping iterations.
    constexpr int bg_min_fit_factor = 3;  // Minimum samples per fitted coefficient factor.

    // Thresholds
    constexpr double source_thresh = 2.0;  // Source-detection SNR threshold.
    constexpr double core_thresh = 4.0;  // Source-core detection threshold.

    // ==========================================
    // Configuration: Stage-3 unbiased noise-stamp quality gates
    // Method: Accept physically compatible candidates with fixed pass/fail cuts before random
    //         selection; never rank candidates by peak, MAD, RMS, or mask fraction.
    // ==========================================
    constexpr double noise_sigma_ratio_min = 0.80;  // Minimum blank-to-source sigma ratio.
    constexpr double noise_sigma_ratio_max = 1.25;  // Maximum blank-to-source sigma ratio.
    constexpr double noise_mad_ratio_min = 0.70;  // Minimum blank-to-source MAD ratio.
    constexpr double noise_mad_ratio_max = 1.30;  // Maximum blank-to-source MAD ratio.
    constexpr double noise_tail_sigma = 2.5;  // Tail-count sigma threshold.
    constexpr double noise_max_tail_fraction = 0.05;  // Maximum blank-stamp tail fraction.
    constexpr double noise_max_mask_fraction = 0.02;  // Maximum blank-stamp masked fraction.

    // ==========================================
    // Configuration: Retained local noise-plane fitting utility
    // Method: Define the square shell and minimum valid fraction exercised by
    //         the standalone NoisePlaneFit regression.
    // ==========================================
    constexpr int noise_region_size = 192;  // Outer local-noise square side length.
    constexpr int noise_inner_size = 96;  // Central exclusion square side length.
    constexpr double noise_plane_min_valid_fraction = 0.30;  // Minimum plane-fit shell fraction.

    // ==========================================
    // Configuration: numerical_fix F6 mode-bar noise-plane estimator
    // Method: Keep these constants identical to f77/sig_para.inc. The estimator obtains robust
    //         block seeds, finds the sky mode and lower-side width, performs two symmetric clipped
    //         plane fits, validates the final plane, then normalizes the amplifier immediately.
    // ==========================================
    constexpr int sig_blocksize = 200;  // Noise-estimator block side length.
    constexpr int sig_block_max = sig_blocksize * sig_blocksize;  // Maximum pixels per block.
    constexpr int sig_max_blocks = 2048;  // Maximum sampled noise blocks.
    constexpr int sig_min_block_pixels = 1000;  // Minimum pixels in one block.
    constexpr int sig_min_block_triples = 1000;  // Minimum valid triples per block.
    constexpr int sig_min_blocks = 4;  // Minimum blocks for a plane fit.
    constexpr int sig_hist_nbin = 256;  // Mode-finding histogram bins.
    constexpr double sig_hist_range = 6.0;  // Histogram range in sigma units.
    constexpr int sig_min_mode_count = 500;  // Minimum samples defining the mode.
    constexpr int sig_min_lower_count = 1000;  // Minimum lower-side samples.
    constexpr double sig_lower_quantile = 0.3173105;  // Lower-side width quantile.
    constexpr double sig_clip_k = 3.0;  // Symmetric clipping sigma.
    constexpr int sig_rdil = 2;  // Pixel stride used by the estimator.
    constexpr int sig_clip_niter = 2;  // Number of clipping iterations.
    constexpr int sig_min_fit_triples = 1000;  // Minimum triples for the final fit.
    constexpr double sig_min_fit_frac = 0.20;  // Minimum retained fit fraction.
    constexpr double sig_median_ratio = 1.2678405;  // Median-to-sigma conversion.
    constexpr double sig_plane_min = 1.0e-8;  // Minimum positive sigma-plane value.
    constexpr double sig_max_plane_ratio = 4.0;  // Maximum plane variation ratio.
    constexpr double sig_pivot_min = 1.0e-8;  // Minimum linear-solve pivot.

    // sig_scale converts the fitted plane into the published 2*sigma^2 convention. Stage two is
    // active; all thresholds and downstream calibrations therefore use an honest sigma convention.
    constexpr double sig_scale_s1 = 0.673475;  // Stage-1 noise calibration candidate.
    constexpr double sig_scale_s2 = 1.027786;  // Stage-2 noise calibration.
    constexpr double sig_scale = sig_scale_s2;  // Active noise calibration selector.

    constexpr int area_max = ns * ns;  // Maximum connected source area.
    constexpr int area_thresh = 6;  // Minimum connected source area.

    constexpr int gal_smooth = 0;  // Galaxy-stamp smoothing radius.
    constexpr int star_smooth = 2;  // Star-stamp smoothing radius.
    constexpr int size_fit_rmax = 4;  // Maximum radius for curvature-size fitting.

    // Point-source Fourier-power morphology statistics. The beta value is a
    // single survey-wide template shape, not a per-source fitted parameter.
    constexpr double point_stat_beta = 0.10;  // Fixed point-source template beta.
    constexpr double point_stat_k_frac = 0.90;  // Fourier-radius fraction used by point statistics.
    constexpr double point_stat_eps = 1.0e-20;  // Numerical denominator floor.
    constexpr double point_stat_min_corr = 1.0e-6;  // Minimum accepted template correlation.

    constexpr double SNR_PSF = 100.0;  // Minimum PSF-star signal-to-noise ratio.

    // Catalogue column indices (shifted to 0-based for C++)
    constexpr int iid = 1 - 1;  // PSF polynomial chi-square field index.
    constexpr int ipixx = 2 - 1;  // Source-center x field index.
    constexpr int ipixy = 3 - 1;  // Source-center y field index.
    constexpr int isig = 4 - 1;  // Local noise sigma field index.
    constexpr int istar = 5 - 1;  // Available PSF-star count field index.
    constexpr int ipeak = 5 - 1;  // Historical peak alias field index.
    constexpr int i_imax = 6 - 1;  // Peak x field index.
    constexpr int i_jmax = 7 - 1;  // Peak y field index.
    constexpr int ih_flux = 8 - 1;  // Half-light flux field index.
    constexpr int ih_area = 9 - 1;  // Source area field index.
    constexpr int iflag = 10 - 1;  // Quality flag field index.
    constexpr int iPSF = 11 - 1;  // Local PSF size field index.
    constexpr int iSNR_F = 12 - 1;  // Fourier SNR field index.
    constexpr int ira = 13 - 1;  // Right-ascension field index.
    constexpr int idec = 14 - 1;  // Declination field index.
    constexpr int igf1 = 15 - 1;  // Field-distortion g1 index.
    constexpr int igf2 = 16 - 1;  // Field-distortion g2 index.
    constexpr int ig1 = 17 - 1;  // Fourier_Quad g1 estimator index.
    constexpr int ig2 = 18 - 1;  // Fourier_Quad g2 estimator index.
    constexpr int ide = 19 - 1;  // Shear response estimator index.
    constexpr int ih1 = 20 - 1;  // Higher-order h1 estimator index.
    constexpr int ih2 = 21 - 1;  // Higher-order h2 estimator index.
    constexpr int icos2 = 22 - 1;  // Spin-2 cosine field index.
    constexpr int isin2 = 23 - 1;  // Spin-2 sine field index.
    constexpr int iparity = 24 - 1;  // WCS parity field index.
    constexpr int igalsizeT = 25 - 1;  // Curvature galaxy-size field index.
    constexpr int ipsfsizeT = 26 - 1;  // Curvature PSF-size field index.
    constexpr int idelta_chi2 = 27 - 1;  // Point-template delta-chi-square field index.
    constexpr int iorth_ext = 28 - 1;  // Orthogonal point-template extension field index.
    constexpr int shear_cat_ncols = iorth_ext + 1;  // Stage-7 shear catalog width.
    constexpr int expo_cat_ncols = shear_cat_ncols + 1;  // Exposure catalog width including chi-square.
    constexpr int ichi2 = shear_cat_ncols;  // Exposure chi-square field index.

    // Band correction parameters
    constexpr double g1_c = 0.0;  // Additive field-distortion g1 correction.
    constexpr double g2_c = 0.0;  // Additive field-distortion g2 correction.
    constexpr double chi2_thresh = 0.1;  // Maximum exposure PSF chi-square.

    // Mathematical constants
    constexpr double pi = 3.14159265358979323846;  // Mathematical pi.
    constexpr double arc_convert = pi / 180.0;  // Degrees-to-radians conversion factor.

    // ==========================================
    static_assert(psf_exposure_min_candidates > 0,
                  "PSF exposure minimum must be positive");
    static_assert(psf_minchi_reference_fraction > 0.0
                      && psf_minchi_reference_fraction <= 1.0,
                  "PSF minChi reference fraction must lie in (0,1]");
    static_assert(psf_minchi_reference_max_per_chip > 0,
                  "PSF minChi reference cap must be positive");
    static_assert(psf_count_pilot_clip_sigma > 0.0,
                  "PSF star-area pilot clipping sigma must be positive");
    static_assert(psf_count_pilot_clip_iterations > 0,
                  "PSF star-area pilot clipping iterations must be positive");
    static_assert(psf_count_zero_mad_quantile >= 0.0
                      && psf_count_zero_mad_quantile < 0.5,
                  "PSF star-area zero-MAD quantile must lie in [0,0.5)");
    static_assert(psf_count_hist_range_sigma > 0.0,
                  "PSF star-area histogram range sigma must be positive");
    static_assert(psf_pair_chi_valid_peak_fraction > 0.0
                      && psf_pair_chi_valid_peak_fraction < 1.0,
                  "PSF pair-chi peak fraction must lie in (0,1)");
    static_assert(psf_bad_fraction_valid_peak_fraction > 0.0
                      && psf_bad_fraction_valid_peak_fraction < 1.0,
                  "PSF bad-pair peak fraction must lie in (0,1)");
    static_assert(psf_type3_elbow_search_height_fraction > 0.0
                      && psf_type3_elbow_search_height_fraction < 1.0,
                  "PSF Type-3 elbow search height fraction must lie in (0,1)");
    static_assert(psf_press_max_removals >= 0,
                  "PSF PRESS removal cap must be non-negative");
    static_assert(psf_loo_min_denom > 0.0 && psf_loo_min_denom < 1.0,
                  "PSF LOO denominator floor must lie in (0,1)");
    static_assert(noise_region_size > noise_inner_size,
                  "noise region must exceed the central exclusion");
    static_assert(noise_inner_size >= nl,
                  "noise inner exclusion must cover the full source extraction region");
    static_assert(noise_region_size % 2 == 0 && noise_inner_size % 2 == 0,
                  "noise region and exclusion sizes must be even");
    static_assert((noise_region_size - noise_inner_size) % 2 == 0,
                  "noise inner exclusion must be centered on the local noise region");
    static_assert(noise_plane_min_valid_fraction > 0.0
                      && noise_plane_min_valid_fraction <= 1.0,
                  "noise plane minimum valid fraction must lie in (0,1]");
}

#endif // LENSING_CONFIG_HPP
