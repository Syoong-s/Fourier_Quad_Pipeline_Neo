#ifndef LENSING_CONFIG_HPP
#define LENSING_CONFIG_HPP

#include <string>
#include <cmath>

namespace LensingConfig {
    // Mathematical constants
    constexpr double pi = 3.14159265358979323846;
    constexpr double arc_convert = pi / 180.0;

    // Image/CCD size parameters
    constexpr int npx = 3000;
    constexpr int npy = 5000;

    // Stage control parameters
    constexpr int ASTROMETRY_trivial = 0;
    constexpr int PROCESS_stage = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23;
    constexpr int include_FLAT = 0;
    constexpr int include_Mask = 2;

    // Catalog/flat paths (originally from para.inc).
    const std::string ASTROMETRY_CAT = "/lustre/home/acct-phyzj/phyzj/jzhang/gaia/gaia_cat_sorted";

    // ==========================================
    // Configuration: Primary external source-catalog directory
    // Method: Seed process_extcat output and process_main input from one mutable path so a
    //         command-line override can update both phases before processing starts.
    // ==========================================
    inline std::string SOURCE_CAT = "/lustre/home/acct-phyzj/share/DES/testy/des_y6_cat";
    const std::string FLAT_PATH = "/lustre/home/acct-phyzj/share/DES/testy/DES_super_flat/i2014";
    const std::string PSF_PATH = "hahahaha";

    // Split parameters
    constexpr int ext_cat = 1;
    constexpr int ext_PSF = 0;
    constexpr int CCD_split = 2;
    constexpr int nct = 12;
    constexpr int ncx = 3;

    // PSF selection and configuration
    constexpr int psf_order = 8;
    constexpr int npo = 64;
    constexpr int npox = 8;
    constexpr int nstar_min = npo * 3 / 2;
   constexpr int npl = 10;
   constexpr int nplx = 2;
   constexpr int nstar_min_local = 16;

    constexpr int step_psf = 100;
    constexpr int deblending = 1;
    constexpr int n_neighbor = 5;

    constexpr int PSF_type = 1;
    constexpr int PSF_Ms = 0;

    // Stamp dimensions
    constexpr int ns = 64;
    constexpr int nsns = ns * ns;
    constexpr int chip_margin = 8;
    constexpr int ns_2 = ns / 2;
    constexpr int nl_2 = ns_2 + chip_margin;
    constexpr int nl = nl_2 * 2;
    constexpr int flag_thresh = 3;
    constexpr int chip_edge_margin = chip_margin;

    constexpr double dz_thresh = 0.1;

    // Catalog sizes and limits
    constexpr int len_g = 40;
    constexpr int len_s = 15;
    constexpr int ngal_max = 4000;
    constexpr int nstar_max = 2000;
    constexpr int npara = 25;
    constexpr int len_sam = 50;

    constexpr int npd = 33;
    constexpr int blocksize = 200;

    // Thresholds
    constexpr double source_thresh = 2.0;
    constexpr double core_thresh = 4.0;
    constexpr double flat_thresh = 0.01;

    // ==========================================
    // Configuration: Function Set_Sig mode-bar noise-plane estimator
    // Method: Keep these constants identical to f77/sig_para.inc. The estimator obtains robust
    //         block seeds, finds the sky mode and lower-side width, performs two symmetric clipped
    //         plane fits, validates the final plane, then normalizes the amplifier immediately.
    // ==========================================
    constexpr int sig_blocksize = 200;
    constexpr int sig_block_max = sig_blocksize * sig_blocksize;
    constexpr int sig_max_blocks = 2048;
    constexpr int sig_min_block_pixels = 1000;
    constexpr int sig_min_block_triples = 1000;
    constexpr int sig_min_blocks = 4;
    constexpr int sig_hist_nbin = 256;
    constexpr double sig_hist_range = 6.0;
    constexpr int sig_min_mode_count = 500;
    constexpr int sig_min_lower_count = 1000;
    constexpr double sig_lower_quantile = 0.3173105;
    constexpr double sig_clip_k = 3.0;
    constexpr int sig_rdil = 2;
    constexpr int sig_clip_niter = 2;
    constexpr int sig_min_fit_triples = 1000;
    constexpr double sig_min_fit_frac = 0.20;
    constexpr double sig_median_ratio = 1.2678405;
    constexpr double sig_plane_min = 1.0e-8;
    constexpr double sig_max_plane_ratio = 4.0;
    constexpr double sig_pivot_min = 1.0e-8;

    // sig_scale converts the fitted plane into the published 2*sigma^2 convention. Stage two is
    // active; all thresholds and downstream calibrations therefore use an honest sigma convention.
    constexpr double sig_scale_s1 = 0.673475;
    constexpr double sig_scale_s2 = 1.027786;
    constexpr double sig_scale = sig_scale_s2;

    constexpr int area_max = ns * ns;
   constexpr int area_thresh = 6;

    constexpr int gal_smooth = 0;
    constexpr int star_smooth = 2;

    constexpr double SNR_PSF = 100.0;
    constexpr double saturation_thresh = 25000.0;

    // Scale conversion
    constexpr double pixel_size = 0.2628; // arcsec

    // Catalogue column indices (shifted to 0-based for C++)
    constexpr int iid = 1 - 1;
    constexpr int ipixx = 2 - 1;
    constexpr int ipixy = 3 - 1;
    constexpr int isig = 4 - 1;
    constexpr int istar = 5 - 1;
    constexpr int ipeak = 5 - 1;
    constexpr int i_imax = 6 - 1;
    constexpr int i_jmax = 7 - 1;
    constexpr int ih_flux = 8 - 1;
    constexpr int ih_area = 9 - 1;
    constexpr int iflag = 10 - 1;
    constexpr int iPSF = 11 - 1;
    constexpr int iSNR_F = 12 - 1;
    constexpr int ira = 13 - 1;
    constexpr int idec = 14 - 1;
    constexpr int igf1 = 15 - 1;
    constexpr int igf2 = 16 - 1;
    constexpr int ig1 = 17 - 1;
    constexpr int ig2 = 18 - 1;
    constexpr int ide = 19 - 1;
    constexpr int ih1 = 20 - 1;
    constexpr int ih2 = 21 - 1;
    constexpr int icos2 = 22 - 1;
    constexpr int isin2 = 23 - 1;
    constexpr int iparity = 24 - 1;
    constexpr int ichi2 = 25 - 1;

    // Max counts
    constexpr int NMAX_EXPO = 25000;
    constexpr int NMAX_CHIP = 62;

    // Band correction parameters
    constexpr double g1_c = -0.001;
    constexpr double g2_c = -0.0003;

    constexpr double chi2_thresh = 0.01;

    // ==================== From cust_para.inc ===============================
    constexpr int chipnx = 2046;
    constexpr int chipny = 4094;

    constexpr double rescale_size = 1.2;

    constexpr int procs_pn = 40;
    constexpr int work_pn = 10;
    constexpr int nblocks = 2;

    constexpr int n_pcs = 100;
    constexpr int npp6th = 28;
    constexpr double pca_negative_eigenvalue_threshold = -1.0e-5;
    constexpr int nmax_star_pchip = 1000000;
}

#endif // LENSING_CONFIG_HPP
