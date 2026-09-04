#ifndef FD_CONFIG_HPP
#define FD_CONFIG_HPP

// ==========================================
// FDConfig — C++ equivalent of the Fortran para.inc
// Method: All FD-test parameters fixed to one set (DES defaults). The two
//         feature switches below select the sigma-estimation method and the
//         star-bar fitting mode at compile time.
// ==========================================

namespace FDConfig {

// ==================== Feature Switches ====================
// Statistical mode: selects how the mean shear (c_best) and its uncertainty
// (sigma) are estimated per spatial bin.
//   PDF_SIGMA          - chi2 sign test + quadratic fitting (c and sigma from PDF)
//   PDF_JACK - PDF chi2 sign test for c_best, jackknife for sigma
//   SWSE_JACK           - SWSE ratio estimator for c, jackknife for sigma
enum class StaticMode {
    PDF_SIGMA,  // PDF chi-square estimate with analytic sigma.
    PDF_JACK,   // PDF chi-square estimate with jackknife sigma.
    SWSE_JACK   // Ratio estimate with jackknife sigma.
};
inline constexpr StaticMode FD_STATIC_MODE = StaticMode::PDF_SIGMA;  // Active FD estimator mode.

// Derived helpers for compile-time branching
inline constexpr bool FD_USE_PDF_STATIS =   // Mode 1 or 2: statis uses PDF chi2 sign test
    (FD_STATIC_MODE == StaticMode::PDF_SIGMA ||
     FD_STATIC_MODE == StaticMode::PDF_JACK);
inline constexpr bool FD_USE_JACKKNIFE =    // Mode 2 or 3: plotComparison does jackknife for sigma
    (FD_STATIC_MODE == StaticMode::PDF_JACK ||
     FD_STATIC_MODE == StaticMode::SWSE_JACK);
inline constexpr bool FD_USE_SWSE_DATA =    // Mode 3: star cut uses SWSE data model
    (FD_STATIC_MODE == StaticMode::SWSE_JACK);

// Star-bar mode: true = per-exposure star bar (NtoN style),
//                false = single global star bar (Nto1 style)
inline constexpr bool FD_PER_EXPOSURE_STAR_BAR = false;  // Fit one star bar per exposure when true.

// ==================== Dimensions ====================
inline constexpr int nmax_per_core = 20000000;  // Maximum sources reserved per MPI rank.
inline constexpr int fd_num = 21;          // spatial bins by field distortion
inline constexpr int PDF_BINS = 4;         // equal-probability inner bins
inline constexpr float gf_lim = 0.0015;      // spatial bin range ±gf_lim
inline constexpr int NMAX = 200;           // fine grid sampling points
inline constexpr int MAX_DUP = 5;          // max duplicate measurements

// ==================== Jackknife / K-means ====================
inline constexpr int N_jack = 50;          // jackknife regions
inline constexpr int nmax_total = 1000000; // max total sources for k-means
inline constexpr int Km_iter = 100;        // k-means iterations

// ==================== Quality-cut thresholds ====================
inline constexpr float snrfcut = 4.0;  // Minimum Fourier signal-to-noise ratio.
inline constexpr float snrlow = 0.0;  // Optional lower source-SNR bound; zero disables it.
inline constexpr float snrhigh = 0.0;  // Optional upper source-SNR bound; zero disables it.
inline constexpr float starcut = 20.0;  // Point-source size cut.
inline constexpr float chi2_thresh = 0.01;  // Maximum exposure chi-square.
inline constexpr float flagcut = 0.0;  // Maximum accepted source quality flag.
inline constexpr float imaxcut = 64.0;  // Maximum source peak x coordinate.
inline constexpr float jmaxcut = 64.0;  // Maximum source peak y coordinate.
inline constexpr float zplow = 0.0;  // Minimum photometric redshift.
inline constexpr float zphigh = 3.0;  // Maximum photometric redshift.
inline constexpr float r_half_thresh = 0.0;  // Optional half-light-radius threshold.
inline constexpr float star_bar_mltp = 3.0;  // Stellar-locus sigma multiplier.
inline constexpr float psf_chi2_mltp = 3.0;  // PSF chi-square sigma multiplier.

// ==================== Star-cut histogram parameters ====================
inline constexpr int n_size_bins = 100;  // Stellar-size histogram bins.
inline constexpr int n_mag_bins = 20;  // Magnitude histogram bins.
inline constexpr float size_min = -2.0;  // Minimum histogram size coordinate.
inline constexpr float size_max = 2.0;  // Maximum histogram size coordinate.
inline constexpr float mag_min_val = 10.0;  // Minimum histogram magnitude.
inline constexpr float mag_max_val = 30.0;  // Maximum histogram magnitude.
inline constexpr int min_bin_count = 100;  // Minimum samples in a usable bin.
inline constexpr float peak_match_tol = 0.05;  // Stellar-peak matching tolerance.
inline constexpr float min_concentration = 0.6;  // Minimum stellar-locus concentration.
inline constexpr float star_phy_min = -0.5;  // Minimum physical stellar size.
inline constexpr float star_phy_max = 0.2;  // Maximum physical stellar size.

// Per-exposure star-cut additional parameters
inline constexpr float stage1_snr = 40.0;  // SNR for histogram accumulation
inline constexpr float stage2_snr = 0.0;  // Second-pass SNR threshold; zero disables it.
inline constexpr float init_win_active = 0.1;  // Active exposure initial size window.
inline constexpr float init_win_fallback = 0.15;  // Fallback exposure initial size window.
inline constexpr float default_s_init = 0.5;  // Default initial stellar-size center.
inline constexpr float clip_nsigma = 3.0;  // Iterative clipping sigma.
inline constexpr float min_clip_limit = 0.015;  // Minimum clipping half-width.
inline constexpr float default_s_std = 0.05;  // Default stellar-size scatter.
inline constexpr float fallback_scut_default = 0.6;  // Fallback stellar-size cut.

// ==================== Bad CCD list (DES) ====================
inline constexpr int bad_ccds[] = {-1};  // DES CCD numbers excluded from analysis.
inline constexpr int n_bad_ccds = 0;  // Number of excluded CCDs.

// ==================== Chip-edge masking ====================
inline constexpr int chip_mask_edge = 50;  // Mask this many pixels from every CCD edge.

}  // namespace FDConfig

#endif  // FD_CONFIG_HPP
