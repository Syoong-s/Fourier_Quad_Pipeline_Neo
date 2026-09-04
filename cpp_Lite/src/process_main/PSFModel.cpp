#include "process_main/PSFModel.hpp"
#include "process_main/ProcessMainState.hpp"
#include "process_main/PSFModelState.hpp"
#include "process_main/PSFCandidateQuality.hpp"
#include "process_main/PSFStarSelection.hpp"
#include "process_main/OutputFile.hpp"
#include "process_main/MPIFailure.hpp"
#include "general/OutputLayout.hpp"
#include "LensingConfig.hpp"
#include "RuntimeConfig.hpp"
#include "process_main/FitsIO.hpp"
#include "process_main/Astrometry.hpp"
#include "general/NumericalRecipes.hpp"
#include "process_main/UniversalUtils.hpp"
#include "process_main/Universalblock.hpp"
#include "process_main/ImageProcessing.hpp"
#include "process_main/ExStar.hpp"
#include "process_main/LinearSolve.hpp"
#include <mpi.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <array>
#include <cstddef>
#include <complex>
#include <initializer_list>
#include <memory>
#include <limits>
#include <new>
#include <stdexcept>

// Extern variables defined elsewhere (e.g. main.cpp)

namespace PSFModel {

    using Internal::ChipPSFState;
    using Internal::CandidatePowerStatus;
    using Internal::ExposurePSFState;

    // Forward declarations of local helper functions
    void readInCandidates(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int& nc, std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state);
    void starSelection(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void applyPressSelection(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void plotStarExpo(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);
    void plotStars(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int nc, const std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state);
    void makePSFLocalFit(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state);

    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics,
        std::vector<double>* leverage);
    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0);

    void getPowerArea(int nx, int ny, const std::vector<float>& power, int& area, float thresh_ratio);
    void getPowerE(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, float thresh_ratio);
    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio);
    void getPSFFWHM(
        const std::vector<float>& power,
        double& FWHM,
        int& star_area);

    // ==========================================
    // Function: Validate one PSF fitting sample
    // Method: Require finite position, shape diagnostics, and every power-spectrum pixel before retaining the star.
    // ==========================================
    static bool isFinitePSFStar(const std::vector<float>& star, std::size_t offset, int pixel_count,
                                double x, double y, double size, double e1, double e2) {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(size) ||
            !std::isfinite(e1) || !std::isfinite(e2) ||
            offset + static_cast<std::size_t>(pixel_count) > star.size()) {
            return false;
        }
        for (int idx = 0; idx < pixel_count; ++idx) {
            if (!std::isfinite(star[offset + idx])) {
                return false;
            }
        }
        return true;
    }

    // ==========================================
    // Function: Escape arbitrary text for an SVG XML text node
    // Method: Replace all five XML-sensitive characters with entity references.
    // ==========================================
    static std::string escapeSVGText(const std::string& text) {
        std::string escaped;
        escaped.reserve(text.size());
        for (const char character : text) {
            switch (character) {
                case '&': escaped += "&amp;"; break;
                case '<': escaped += "&lt;"; break;
                case '>': escaped += "&gt;"; break;
                case '"': escaped += "&quot;"; break;
                case '\'': escaped += "&apos;"; break;
                default: escaped += character; break;
            }
        }
        return escaped;
    }

    // ==========================================
    // Function: Write one exposure-level integer-area locus SVG diagnostic
    // Method: Render the same-pass two-count histograms, pilot/peak/MAD/elbow
    //         decisions, and final guarded cuts before PRESS rejection.
    // ==========================================
    static void writePSFCountLocusSVG(
        const std::string& dirOutput,
        const std::string& exposure,
        const Internal::PSFCountLocus& locus,
        const Internal::PSFCountLocusDiagnostics& diagnostics) {
        if (!locus.valid || diagnostics.histogram.empty()
            || diagnostics.histogram.size()
                != diagnostics.smoothed_histogram.size()
            || diagnostics.peak_bin < 0
            || diagnostics.peak_bin
                >= static_cast<int>(diagnostics.histogram.size())) {
            return;
        }
        const std::size_t bin_count = diagnostics.histogram.size();
        const int first_count = diagnostics.histogram_first_count;
        const int last_count = diagnostics.histogram_last_count;
        const int expected_bin_count = last_count >= first_count
            ? (last_count - first_count)
                / Internal::PSFCountHistogramBinWidth + 1
            : 0;
        if (expected_bin_count != static_cast<int>(bin_count)
            || !std::isfinite(diagnostics.mad_lower)
            || !std::isfinite(diagnostics.mad_upper)) {
            return;
        }
        const double histogram_lower =
            static_cast<double>(first_count) - 0.5;
        const double histogram_upper =
            static_cast<double>(first_count)
                + static_cast<double>(Internal::PSFCountHistogramBinWidth)
                    * static_cast<double>(bin_count) - 0.5;
        const bool has_left_elbow = diagnostics.left_elbow_bin >= 0
            && diagnostics.left_elbow_bin < static_cast<int>(bin_count);
        const bool has_right_elbow = diagnostics.right_elbow_bin >= 0
            && diagnostics.right_elbow_bin < static_cast<int>(bin_count);
        const double left_elbow_value = has_left_elbow
            ? Internal::psfCountHistogramBinCenter(
                first_count, diagnostics.left_elbow_bin)
            : 0.0;
        const double right_elbow_value = has_right_elbow
            ? Internal::psfCountHistogramBinCenter(
                first_count, diagnostics.right_elbow_bin)
            : 0.0;
        const bool has_gaia_histogram =
            !diagnostics.gaia_histogram.empty()
            && diagnostics.gaia_histogram.size()
                == diagnostics.histogram.size();
        const bool has_minchi_survivor_histogram =
            !diagnostics.minchi_survivor_histogram.empty()
            && diagnostics.minchi_survivor_histogram.size()
                == diagnostics.histogram.size();
        const bool has_selected_group_histogram =
            !diagnostics.selected_group_histogram.empty()
            && diagnostics.selected_group_histogram.size()
                == diagnostics.histogram.size();
        std::ostringstream quantile_range_mode;
        quantile_range_mode << "Q("
                            << LensingConfig::psf_count_zero_mad_quantile
                            << ")-Q("
                            << 1.0 - LensingConfig::psf_count_zero_mad_quantile
                            << ")";

        constexpr double canvas_width = 1200.0;
        constexpr double canvas_height = 1160.0;
        constexpr double plot_left = 90.0;
        constexpr double plot_right = 870.0;
        constexpr double plot_top = 110.0;
        constexpr double plot_bottom = 650.0;
        const double plot_width = plot_right - plot_left;
        const double plot_height = plot_bottom - plot_top;

        double x_min = std::min(diagnostics.pilot_lower, locus.lower);
        double x_max = std::max(diagnostics.pilot_upper, locus.upper);
        x_min = std::min(x_min, histogram_lower);
        x_max = std::max(x_max, histogram_upper);
        x_min = std::min(x_min, diagnostics.pilot_center);
        x_max = std::max(x_max, diagnostics.pilot_center);
        x_min = std::min(x_min, diagnostics.mad_lower);
        x_max = std::max(x_max, diagnostics.mad_upper);
        if (has_left_elbow) x_min = std::min(x_min, left_elbow_value);
        if (has_right_elbow) x_max = std::max(x_max, right_elbow_value);
        if (diagnostics.has_gaia_median) {
            x_min = std::min(x_min, diagnostics.gaia_median);
            x_max = std::max(x_max, diagnostics.gaia_median);
        }
        const double x_padding = std::max((x_max - x_min) * 0.05, 1.0e-9);
        x_min -= x_padding;
        x_max += x_padding;

        double y_max = 0.0;
        for (const double count : diagnostics.histogram) {
            y_max = std::max(y_max, count);
        }
        for (const double count : diagnostics.smoothed_histogram) {
            y_max = std::max(y_max, count);
        }
        if (has_gaia_histogram) {
            for (const double count : diagnostics.gaia_histogram) {
                y_max = std::max(y_max, count);
            }
        }
        if (has_minchi_survivor_histogram) {
            for (const double count : diagnostics.minchi_survivor_histogram) {
                y_max = std::max(y_max, count);
            }
        }
        if (has_selected_group_histogram) {
            for (const double count : diagnostics.selected_group_histogram) {
                y_max = std::max(y_max, count);
            }
        }
        y_max = std::max(1.0, y_max * 1.08);

        const auto mapX = [&](double value) {
            return plot_left + (value - x_min) / (x_max - x_min) * plot_width;
        };
        const auto mapY = [&](double value) {
            return plot_bottom - value / y_max * plot_height;
        };

        const std::string filename = dirOutput + "/stamps/svg_StarLocus/"
            + exposure + "_locus.svg";
        MainIO::OutputFile output(filename);
        output << std::fixed << std::setprecision(6);
        output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
               << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
               << canvas_width << "\" height=\"" << canvas_height
               << "\" viewBox=\"0 0 " << canvas_width << ' ' << canvas_height
               << "\">\n"
               << "  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n"
               << "  <text x=\"600\" y=\"38\" text-anchor=\"middle\" "
                  "font-family=\"sans-serif\" font-size=\"24\" font-weight=\"bold\">"
                  "PSF Star Pixel-Count Locus</text>\n"
               << "  <text x=\"600\" y=\"68\" text-anchor=\"middle\" "
                  "font-family=\"sans-serif\" font-size=\"16\">Exposure: "
               << escapeSVGText(exposure) << "</text>\n";

        for (int tick = 0; tick <= 5; ++tick) {
            const double fraction = static_cast<double>(tick) / 5.0;
            const double x_value = x_min + fraction * (x_max - x_min);
            const double x = mapX(x_value);
            const double y_value = fraction * y_max;
            const double y = mapY(y_value);
            output << "  <line x1=\"" << x << "\" y1=\"" << plot_top
                   << "\" x2=\"" << x << "\" y2=\"" << plot_bottom
                   << "\" stroke=\"#eeeeee\"/>\n"
                   << "  <text x=\"" << x << "\" y=\"" << plot_bottom + 25.0
                   << "\" text-anchor=\"middle\" font-family=\"sans-serif\" "
                      "font-size=\"12\">" << x_value << "</text>\n"
                   << "  <line x1=\"" << plot_left << "\" y1=\"" << y
                   << "\" x2=\"" << plot_right << "\" y2=\"" << y
                   << "\" stroke=\"#eeeeee\"/>\n"
                   << "  <text x=\"" << plot_left - 12.0 << "\" y=\"" << y + 4.0
                   << "\" text-anchor=\"end\" font-family=\"sans-serif\" "
                      "font-size=\"12\">" << y_value << "</text>\n";
        }

        output << "  <g id=\"raw-histogram\" fill=\"#b8bec7\">\n";
        for (std::size_t bin = 0; bin < bin_count; ++bin) {
            const double left_value = static_cast<double>(first_count)
                + static_cast<double>(Internal::PSFCountHistogramBinWidth)
                    * static_cast<double>(bin) - 0.5;
            const double right_value = left_value
                + static_cast<double>(Internal::PSFCountHistogramBinWidth);
            const double left = mapX(left_value);
            const double right = mapX(right_value);
            const double top = mapY(diagnostics.histogram[bin]);
            output << "    <rect x=\"" << left << "\" y=\"" << top
                   << "\" width=\"" << std::max(0.0, right - left - 0.5)
                   << "\" height=\"" << plot_bottom - top << "\"/>\n";
        }
        output << "  </g>\n"
               << "  <polyline id=\"smoothed-histogram\" fill=\"none\" "
                  "stroke=\"#2468b4\" stroke-width=\"3\" points=\"";
        for (std::size_t bin = 0; bin < bin_count; ++bin) {
            const double center = Internal::psfCountHistogramBinCenter(
                first_count, static_cast<int>(bin));
            output << mapX(center) << ','
                   << mapY(diagnostics.smoothed_histogram[bin]) << ' ';
        }
        output << "\"/>\n";
        if (has_gaia_histogram) {
            output << "  <polyline id=\"gaia-histogram\" fill=\"none\" "
                      "stroke=\"#2ca02c\" stroke-width=\"2\" points=\"";
            for (std::size_t bin = 0; bin < bin_count; ++bin) {
                const double center = Internal::psfCountHistogramBinCenter(
                    first_count, static_cast<int>(bin));
                output << mapX(center) << ','
                       << mapY(diagnostics.gaia_histogram[bin]) << ' ';
            }
            output << "\"/>\n";
        }
        if (has_minchi_survivor_histogram) {
            output << "  <polyline id=\"minchi-survivor-histogram\" fill=\"none\" "
                      "stroke=\"#e377c2\" stroke-width=\"3\" points=\"";
            for (std::size_t bin = 0; bin < bin_count; ++bin) {
                const double center = Internal::psfCountHistogramBinCenter(
                    first_count, static_cast<int>(bin));
                output << mapX(center) << ','
                       << mapY(diagnostics.minchi_survivor_histogram[bin]) << ' ';
            }
            output << "\"/>\n";
        }
        if (has_selected_group_histogram) {
            output << "  <polyline id=\"selected-group-histogram\" fill=\"none\" "
                      "stroke=\"#17becf\" stroke-width=\"3\" points=\"";
            for (std::size_t bin = 0; bin < bin_count; ++bin) {
                const double center = Internal::psfCountHistogramBinCenter(
                    first_count, static_cast<int>(bin));
                output << mapX(center) << ','
                       << mapY(diagnostics.selected_group_histogram[bin]) << ' ';
            }
            output << "\"/>\n";
        }

        const double peak_value = Internal::psfCountHistogramBinCenter(
            first_count, diagnostics.peak_bin);
        output << "  <line id=\"selected-peak\" x1=\"" << mapX(peak_value)
               << "\" y1=\"" << plot_top << "\" x2=\"" << mapX(peak_value)
               << "\" y2=\"" << plot_bottom
               << "\" stroke=\"#f28e2b\" stroke-width=\"2\" "
                  "stroke-dasharray=\"3,3\"/>\n";
        output << "  <line id=\"pilot-center\" x1=\""
               << mapX(diagnostics.pilot_center) << "\" y1=\"" << plot_top
               << "\" x2=\"" << mapX(diagnostics.pilot_center)
               << "\" y2=\"" << plot_bottom
               << "\" stroke=\"#9467bd\" stroke-width=\"2\" "
                  "stroke-dasharray=\"6,3\"/>\n";
        if (diagnostics.has_gaia_median) {
            output << "  <line id=\"gaia-median\" x1=\""
                   << mapX(diagnostics.gaia_median) << "\" y1=\"" << plot_top
                   << "\" x2=\"" << mapX(diagnostics.gaia_median)
                   << "\" y2=\"" << plot_bottom
                   << "\" stroke=\"#2ca02c\" stroke-width=\"2\" "
                      "stroke-dasharray=\"8,4,2,4\"/>\n";
        }
        output << "  <line id=\"mad-lower\" x1=\""
               << mapX(diagnostics.mad_lower)
               << "\" y1=\"" << plot_top << "\" x2=\""
               << mapX(diagnostics.mad_lower) << "\" y2=\"" << plot_bottom
               << "\" stroke=\"#8c564b\" stroke-width=\"2\" "
                  "stroke-dasharray=\"5,3\"/>\n"
               << "  <line id=\"mad-upper\" x1=\""
               << mapX(diagnostics.mad_upper)
               << "\" y1=\"" << plot_top << "\" x2=\""
               << mapX(diagnostics.mad_upper) << "\" y2=\"" << plot_bottom
               << "\" stroke=\"#8c564b\" stroke-width=\"2\" "
                  "stroke-dasharray=\"5,3\"/>\n";
        if (has_left_elbow) {
            output << "  <line id=\"left-elbow\" x1=\""
                   << mapX(left_elbow_value) << "\" y1=\"" << plot_top
                   << "\" x2=\"" << mapX(left_elbow_value) << "\" y2=\""
                   << plot_bottom << "\" stroke=\"#bc8f00\" stroke-width=\"2\" "
                      "stroke-dasharray=\"2,3\"/>\n";
        }
        if (has_right_elbow) {
            output << "  <line id=\"right-elbow\" x1=\""
                   << mapX(right_elbow_value) << "\" y1=\"" << plot_top
                   << "\" x2=\"" << mapX(right_elbow_value) << "\" y2=\""
                   << plot_bottom << "\" stroke=\"#bc8f00\" stroke-width=\"2\" "
                      "stroke-dasharray=\"2,3\"/>\n";
        }
        output << "  <line id=\"locus-lower\" x1=\"" << mapX(locus.lower)
               << "\" y1=\"" << plot_top << "\" x2=\"" << mapX(locus.lower)
               << "\" y2=\"" << plot_bottom
               << "\" stroke=\"#d62728\" stroke-width=\"3\" "
                  "stroke-dasharray=\"8,5\"/>\n"
               << "  <line id=\"locus-upper\" x1=\"" << mapX(locus.upper)
               << "\" y1=\"" << plot_top << "\" x2=\"" << mapX(locus.upper)
               << "\" y2=\"" << plot_bottom
               << "\" stroke=\"#d62728\" stroke-width=\"3\" "
                  "stroke-dasharray=\"8,5\"/>\n"
               << "  <line id=\"locus-center\" x1=\"" << mapX(locus.center)
               << "\" y1=\"" << plot_top << "\" x2=\"" << mapX(locus.center)
               << "\" y2=\"" << plot_bottom
               << "\" stroke=\"#111111\" stroke-width=\"3\"/>\n"
               << "  <line x1=\"" << plot_left << "\" y1=\"" << plot_bottom
               << "\" x2=\"" << plot_right << "\" y2=\"" << plot_bottom
               << "\" stroke=\"black\" stroke-width=\"2\"/>\n"
               << "  <line x1=\"" << plot_left << "\" y1=\"" << plot_top
               << "\" x2=\"" << plot_left << "\" y2=\"" << plot_bottom
               << "\" stroke=\"black\" stroke-width=\"2\"/>\n"
               << "  <text x=\"" << (plot_left + plot_right) / 2.0
               << "\" y=\"715\" text-anchor=\"middle\" font-family=\"sans-serif\" "
                  "font-size=\"16\">exp(-1) Pixel Count</text>\n"
               << "  <text x=\"25\" y=\"380\" text-anchor=\"middle\" "
                  "transform=\"rotate(-90 25 380)\" font-family=\"sans-serif\" "
                  "font-size=\"16\">Candidate count</text>\n";

        output << std::setprecision(10);
        output << "  <g font-family=\"sans-serif\" font-size=\"14\" fill=\"#222222\">\n"
               << "    <text x=\"900\" y=\"115\">Samples = "
               << diagnostics.sample_count << "</text>\n"
               << "    <text x=\"900\" y=\"138\">Pilot source = "
               << (diagnostics.pilot_uses_gaia ? "Gaia" : "all")
               << "</text>\n"
               << "    <text x=\"900\" y=\"161\">Pilot samples = "
               << diagnostics.pilot_retained_count << " / "
               << diagnostics.pilot_input_count << "</text>\n"
               << "    <text x=\"900\" y=\"184\">Pilot center = "
               << diagnostics.pilot_center << "</text>\n"
               << "    <text x=\"900\" y=\"207\">Pilot width = "
               << diagnostics.pilot_width << "</text>\n"
               << "    <text x=\"900\" y=\"230\">Pilot lower = "
               << diagnostics.pilot_lower << "</text>\n"
               << "    <text x=\"900\" y=\"253\">Pilot upper = "
               << diagnostics.pilot_upper << "</text>\n"
               << "    <text x=\"900\" y=\"276\">Pilot range mode = "
               << (diagnostics.pilot_uses_quantile_range
                       ? quantile_range_mode.str() : "MAD")
               << "</text>\n"
               << "    <text x=\"900\" y=\"299\">Zero-MAD clip rejected = "
               << (diagnostics.pilot_rejected_zero_mad_clip ? "yes" : "no")
               << "</text>\n"
               << "    <text x=\"900\" y=\"322\">Histogram samples = "
               << diagnostics.histogram_sample_count << "</text>\n"
               << "    <text x=\"900\" y=\"345\">Below / above = "
               << diagnostics.histogram_below_count << " / "
               << diagnostics.histogram_above_count << "</text>\n"
               << "    <text x=\"900\" y=\"376\">Final center = "
               << locus.center << "</text>\n"
               << "    <text x=\"900\" y=\"399\">Lower width = "
               << locus.lower_width << "</text>\n"
               << "    <text x=\"900\" y=\"422\">Upper width = "
               << locus.upper_width << "</text>\n"
               << "    <text x=\"900\" y=\"445\">Final sigma = "
               << LensingConfig::psf_count_locus_sigma << "</text>\n"
               << "    <text x=\"900\" y=\"468\">MAD lower = "
               << diagnostics.mad_lower << "</text>\n"
               << "    <text x=\"900\" y=\"491\">MAD upper = "
               << diagnostics.mad_upper << "</text>\n"
               << "    <text x=\"900\" y=\"514\">Final lower = "
               << locus.lower << "</text>\n"
               << "    <text x=\"900\" y=\"537\">Final upper = "
               << locus.upper << "</text>\n"
               << "    <text x=\"900\" y=\"560\">Gaia matches = "
               << diagnostics.gaia_match_count << "</text>\n"
               << "    <text x=\"900\" y=\"583\">Gaia histogram = "
               << diagnostics.gaia_histogram_sample_count << "</text>\n"
               << "    <text x=\"900\" y=\"606\">Gaia below / above = "
               << diagnostics.gaia_histogram_below_count << " / "
               << diagnostics.gaia_histogram_above_count << "</text>\n"
               << "    <text x=\"900\" y=\"629\">MinChi survivors = "
               << diagnostics.minchi_survivor_count << "</text>\n"
               << "    <text x=\"900\" y=\"652\">Pre-PRESS selected = "
               << diagnostics.selected_group_count << "</text>\n";
        if (diagnostics.has_gaia_median) {
            output << "    <text x=\"900\" y=\"675\">Gaia raw median = "
                   << diagnostics.gaia_median << "</text>\n";
        }
        output << "    <text x=\"900\" y=\"698\">Left elbow = ";
        if (has_left_elbow) output << left_elbow_value;
        else output << "unavailable";
        output << "</text>\n"
               << "    <text x=\"900\" y=\"721\">Right elbow = ";
        if (has_right_elbow) output << right_elbow_value;
        else output << "unavailable";
        output << "</text>\n"
               << "    <text x=\"900\" y=\"744\">Left guard applied = "
               << (diagnostics.left_elbow_guard_applied ? "yes" : "no")
               << "</text>\n"
               << "    <text x=\"900\" y=\"767\">Right guard applied = "
               << (diagnostics.right_elbow_guard_applied ? "yes" : "no")
               << "</text>\n"
               << "    <text x=\"930\" y=\"820\">raw histogram</text>\n"
               << "    <text x=\"930\" y=\"846\">smoothed histogram</text>\n";
        if (has_gaia_histogram) {
            output << "    <text x=\"930\" y=\"872\">Gaia histogram</text>\n";
        }
        if (has_minchi_survivor_histogram) {
            output << "    <text x=\"930\" y=\"898\">MinChi survivors</text>\n";
        }
        if (has_selected_group_histogram) {
            output << "    <text x=\"930\" y=\"924\">Pre-PRESS selected</text>\n";
        }
        output << "    <text x=\"930\" y=\"950\">selected peak</text>\n"
               << "    <text x=\"930\" y=\"976\">pilot center</text>\n"
               << "    <text x=\"930\" y=\"1002\">locus center</text>\n"
               << "    <text x=\"930\" y=\"1028\">pre-guard MAD cuts</text>\n"
               << "    <text x=\"930\" y=\"1054\">outer elbows</text>\n"
               << "    <text x=\"930\" y=\"1080\">final guarded cuts</text>\n";
        if (diagnostics.has_gaia_median) {
            output << "    <text x=\"930\" y=\"1106\">Gaia raw median</text>\n";
        }
        output << "  </g>\n"
               << "  <rect x=\"900\" y=\"807\" width=\"20\" height=\"14\" "
                  "fill=\"#b8bec7\"/>\n"
               << "  <line x1=\"900\" y1=\"841\" x2=\"920\" y2=\"841\" "
                  "stroke=\"#2468b4\" stroke-width=\"3\"/>\n";
        if (has_gaia_histogram) {
            output << "  <line x1=\"900\" y1=\"867\" x2=\"920\" y2=\"867\" "
                      "stroke=\"#2ca02c\" stroke-width=\"2\"/>\n";
        }
        if (has_minchi_survivor_histogram) {
            output << "  <line x1=\"900\" y1=\"893\" x2=\"920\" y2=\"893\" "
                      "stroke=\"#e377c2\" stroke-width=\"3\"/>\n";
        }
        if (has_selected_group_histogram) {
            output << "  <line x1=\"900\" y1=\"919\" x2=\"920\" y2=\"919\" "
                      "stroke=\"#17becf\" stroke-width=\"3\"/>\n";
        }
        output << "  <line x1=\"900\" y1=\"945\" x2=\"920\" y2=\"945\" "
                  "stroke=\"#f28e2b\" stroke-width=\"2\" stroke-dasharray=\"3,3\"/>\n"
               << "  <line x1=\"900\" y1=\"971\" x2=\"920\" y2=\"971\" "
                  "stroke=\"#9467bd\" stroke-width=\"2\" stroke-dasharray=\"6,3\"/>\n"
               << "  <line x1=\"900\" y1=\"997\" x2=\"920\" y2=\"997\" "
                  "stroke=\"#111111\" stroke-width=\"3\"/>\n"
               << "  <line x1=\"900\" y1=\"1023\" x2=\"920\" y2=\"1023\" "
                  "stroke=\"#8c564b\" stroke-width=\"2\" stroke-dasharray=\"5,3\"/>\n"
               << "  <line x1=\"900\" y1=\"1049\" x2=\"920\" y2=\"1049\" "
                  "stroke=\"#bc8f00\" stroke-width=\"2\" stroke-dasharray=\"2,3\"/>\n"
               << "  <line x1=\"900\" y1=\"1075\" x2=\"920\" y2=\"1075\" "
                  "stroke=\"#d62728\" stroke-width=\"3\" stroke-dasharray=\"8,5\"/>\n";
        if (diagnostics.has_gaia_median) {
            output << "  <line x1=\"900\" y1=\"1101\" x2=\"920\" y2=\"1101\" "
                      "stroke=\"#2ca02c\" stroke-width=\"2\" "
                      "stroke-dasharray=\"8,4,2,4\"/>\n";
        }
        output << "</svg>\n";
    }

    // ==========================================
    // Function: Compute a checked runtime allocation size
    // Method: Multiply dimensions in size_t and abort coherently before an
    //         overflow could under-allocate a shared PSF buffer.
    // ==========================================
    static std::size_t checkedElementCount(
        std::initializer_list<std::size_t> dimensions,
        const std::string& operation) {
        std::size_t count = 1;
        for (const std::size_t dimension : dimensions) {
            if (dimension != 0 &&
                count > std::numeric_limits<std::size_t>::max() / dimension) {
                MPIFailure::abortWorld(operation, "runtime dimensions overflow size_t");
            }
            count *= dimension;
        }
        return count;
    }

    // ==========================================
    // Function: Read one PSF stamp cube with its catalog-defined shape
    // Method: Decode the FITS axes first and abort the MPI world unless width,
    //         height, and plane count exactly match the direct consumer.
    // ==========================================
    static void readExpectedStampCube(const std::string& filename, int nx, int ny,
                                      int count, std::vector<float>& stamps,
                                      const std::string& operation) {
        FitsIO::StampCubeShape shape;
        if (FitsIO::readStampCube(filename, shape, stamps)
            && shape.matches(nx, ny, count)) {
            return;
        }

        std::ostringstream detail;
        detail << filename << " expected=" << nx << "x" << ny << "x" << count
               << " actual=" << shape.nx << "x" << shape.ny << "x"
               << shape.count;
        MPIFailure::abortWorld(operation, detail.str());
    }

    // ==========================================
    // Function: Run Stage-5 Lite PSF modeling for one exposure
    // Method: Execute Gaia/local candidate selection, analytic PRESS, diagnostics,
    //         and the sole normalized polynomial fit with runtime geometry.
    // ==========================================
    void procPSF(int iexpo) {
        if (iexpo <= 0 || iexpo > static_cast<int>(ProcessMain::state.exposure_files.size())) {
            std::cerr << "Error: invalid iexpo index: " << iexpo << std::endl;
            return;
        }
        std::string expo_file_path = ProcessMain::state.exposure_files[iexpo - 1];
        std::vector<std::string> imageFiles;
        std::string dirOutput;
        UniversalUtils::getImageList(expo_file_path, imageFiles, dirOutput);

        const int nchip = static_cast<int>(imageFiles.size());
        const LensingRuntimeConfig& lensing =
            RuntimeConfigStore::get().lensing;
        if (nchip > lensing.nmax_chip) {
            MPIFailure::abortWorld(
                "run Lite Stage-5 PSF model",
                "exposure chip count " + std::to_string(nchip)
                    + " exceeds configured nmax_chip="
                    + std::to_string(lensing.nmax_chip));
        }

        auto state_ptr = std::make_unique<ExposurePSFState>(nchip);
        ExposurePSFState& state = *state_ptr;

        int nc = 0;
        std::vector<std::array<double, 4>> p_chip(
            static_cast<std::size_t>(lensing.nmax_chip),
            {0.0, 0.0, 0.0, 0.0});

        readInCandidates(nchip, imageFiles, dirOutput, nc, p_chip, state);
        starSelection(nchip, imageFiles, dirOutput, state);
        applyPressSelection(nchip, imageFiles, dirOutput, state);
        plotStarExpo(nchip, imageFiles, dirOutput, state);
        plotStars(nchip, imageFiles, dirOutput, nc, p_chip, state);
        makePSFLocalFit(nchip, imageFiles, dirOutput, state);
    }

    // Local Helper Routines

    // ==========================================
    // Function: Load star candidates and Fourier-power stamps for one exposure
    // Method: Preserve runtime cube I/O and chip geometry while rejecting
    //         invalid corrected spectra before diagnostics or pairwise chi.
    // ==========================================
    void readInCandidates(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int& nc, std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state) {
        int ns = LensingConfig::ns;
        const LensingRuntimeConfig& lensing =
            RuntimeConfigStore::get().lensing;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string headname = dirOutput + "/astrometry/Head/" + prefix_e + ".head";

        nc = 0;

        for (int k = 0; k < nchip; ++k) {
            ChipPSFState& chip = state.chips[k];
            chip.stars.clear();
            chip.stars.reserve(LensingConfig::nstar_max);
            chip.selection.clear();
            chip.fit.clear();

            const Universalblock::NormStatus normStatus =
                Universalblock::checkNorm(imageFiles[k], dirOutput);
            if (normStatus == Universalblock::NormStatus::Invalid) {
                continue;
            }
            if (normStatus != Universalblock::NormStatus::Valid) {
                MPIFailure::abortWorld(
                    "validate PSF chip norm",
                    Universalblock::normErrorDetail(
                        normStatus, imageFiles[k], dirOutput));
            }

            double cRPIX[2] = {0.0, 0.0};
            double cD[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
            double cRVAL[2] = {0.0, 0.0};
            double PU[2][LensingConfig::npd];
            int ierror = 0;

            Astrometry::readAstrometryPara(headname, k + 1, cRPIX, cD, cRVAL, PU, LensingConfig::npd, ierror);

            if (ierror == 1) continue;

            nc++;
            double x = 1.0;
            double y = 1.0;
            double xx = 0.0, yy = 0.0;
            Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
            p_chip[nc - 1][0] = xx;
            p_chip[nc - 1][1] = yy;

            x = static_cast<double>(lensing.chipnx);
            y = static_cast<double>(lensing.chipny);
            Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
            p_chip[nc - 1][2] = xx;
            p_chip[nc - 1][3] = yy;

            std::string prefix = UniversalUtils::getPrefix(imageFiles[k]);
            std::string filepath = OutputLayout::chipPath(
                dirOutput, "stamps/dat_StarCanInfo", prefix, "_star_can_info.dat");

            std::ifstream infile(filepath);
            if (!infile.is_open()) {
                MPIFailure::abortWorld("read PSF star-candidate info", filepath);
            }

            std::string header;
            if (!std::getline(infile, header)) {
                MPIFailure::abortWorld(
                    "read PSF star-candidate header", filepath);
            }

            std::string line;
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                float aa[4] = {0.0f, 0.0f, 0.0f, 0.0f};
                if (iss >> aa[0] >> aa[1] >> aa[2] >> aa[3]) {
                    ChipPSFState::StarRow row{};
                    row[0] = aa[0];
                    row[1] = aa[1];
                    row[2] = aa[2];
                    row[3] = aa[3];
                    chip.stars.push_back(row);
                }
            }
            infile.close();

            const int nstar = state.getNStar(k);
            chip.selection.assign(
                static_cast<std::size_t>(nstar), Internal::StarSelectionState{});

            if (nstar > 0) {
                std::string stampPath = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");
                std::vector<float> star;
                readExpectedStampCube(
                    stampPath, ns, ns, nstar, star,
                    "read PSF star-candidate power cube");

                for (int i = 0; i < nstar; ++i) {
                    state.getStarPara(k, i, 4) = 1.0;

                    std::vector<float> source_p(ns * ns);
                    for (int v = 0; v < ns; ++v) {
                        for (int u = 0; u < ns; ++u) {
                            source_p[v * ns + u] = star[static_cast<size_t>(i) * ns * ns + v * ns + u];
                        }
                    }

                    double sum_power = 0.0;
                    double chi_window_sum = 0.0;
                    if (Internal::assessCandidatePower(
                            ns, ns, source_p, sum_power, chi_window_sum)
                        != CandidatePowerStatus::Accepted) {
                        state.getStarPara(k, i, 4) = -1.0;
                        continue;
                    }

                    x = state.getStarPara(k, i, 1);
                    y = state.getStarPara(k, i, 2);
                    Astrometry::xyToXxyy(x, y, xx, yy, cRPIX, cD);
                    state.getStarPara(k, i, 5) = xx;
                    state.getStarPara(k, i, 6) = yy;

                    std::array<double, 2> ee = {0.0, 0.0};
                    double size = 0.0;
                    getPowerAll(ns, ns, source_p, ee, size, 0.02f);
                    state.getStarPara(k, i, 7) = size;
                    state.getStarPara(k, i, 8) = ee[0];
                    state.getStarPara(k, i, 9) = ee[1];

                    double FWHM = 0.0;
                    int star_area = 0;
                    getPSFFWHM(source_p, FWHM, star_area);
                    if (!Internal::candidateDiagnosticsAreFinite(
                            size, ee[0], ee[1])
                        || star_area <= 0) {
                        state.getStarPara(k, i, 4) = -1.0;
                        continue;
                    }
                    state.getStarPara(k, i, 10) = FWHM;
                    state.getStarPara(k, i, 11) = 1.0 / sum_power;
                    state.getStarPara(
                        k, i, ChipPSFState::star_area_index) =
                            static_cast<double>(star_area);

                    Internal::StarSelectionState& selection = chip.selection[i];
                    selection.full_power_sum = sum_power;
                    selection.chi_window_sum = chi_window_sum;
                    const Internal::PSFChiWindow window =
                        Internal::getPSFChiWindow(ns);
                    selection.chi_window.reserve(
                        static_cast<std::size_t>(window.pixelCount()));
                    for (int row = window.first; row <= window.last; ++row) {
                        for (int column = window.first;
                             column <= window.last; ++column) {
                            selection.chi_window.push_back(
                                source_p[row * ns + column]);
                        }
                    }
                }


            }
        }
    }

    // ==========================================
    // Function: Reject every PSF candidate in an exposure
    // Method: Clear all scientific selection flags and release cached windows.
    // ==========================================
    static void rejectExposureCandidates(ExposurePSFState& state) {
        for (ChipPSFState& chip : state.chips) {
            chip.fit.clear();
            for (std::size_t index = 0; index < chip.stars.size(); ++index) {
                chip.stars[index][4] = -1.0;
                if (index >= chip.selection.size()) continue;
                Internal::StarSelectionState& selection = chip.selection[index];
                selection.in_size_locus = false;
                selection.selected_group = false;
                selection.selected_press = false;
                selection.bad_pair_fraction = 0.0;
                std::vector<float>().swap(selection.chi_window);
            }
        }
    }

    // ==========================================
    // Function: Estimate a project-compatible upper-tail threshold
    // Method: Use the legacy low-side peak/width estimator when finite, fall
    //         back to robust quantiles, and keep all values when width vanishes.
    // ==========================================
    static float estimateUpperTailThreshold(
        const std::vector<float>& input,
        double sigma_cut) {
        std::vector<float> values;
        values.reserve(input.size());
        for (float value : input) {
            if (std::isfinite(value)) values.push_back(value);
        }
        if (values.size() <= 4 || sigma_cut <= 0.0) {
            return std::numeric_limits<float>::infinity();
        }

        float peak = 0.0f;
        float width = 0.0f;
        NumericalRecipes::getPeakWidthLowSide(values, peak, width);
        if (std::isfinite(peak) && std::isfinite(width) && width > 0.0f) {
            return peak + static_cast<float>(sigma_cut) * width;
        }

        UniversalUtils::getMedSig(values, peak, width);
        if (std::isfinite(peak) && std::isfinite(width) && width > 0.0f) {
            return peak + static_cast<float>(sigma_cut) * width;
        }
        return *std::max_element(values.begin(), values.end());
    }

    // ==========================================
    // Function: Read one chip's matched Gaia image positions
    // Method: Open the existing astro product once, accept zero matches, and
    //         abort on missing or malformed Gaia-mode Stage-2 data.
    // ==========================================
    static std::vector<std::array<double, 2>> readAstrometryGaiaPositions(
        const std::string& imageFile,
        const std::string& dirOutput) {
        std::vector<std::array<double, 2>> gaia_xy;
        const std::string prefix = UniversalUtils::getPrefix(imageFile);
        const std::string filename = OutputLayout::chipPath(
            dirOutput, "astrometry/dat_Astro", prefix, "_astro.dat");
        std::ifstream input(filename);
        if (!input.is_open()) {
            MPIFailure::abortWorld("read PSF Gaia astrometry matches", filename);
        }
        std::string error;
        const Internal::AstrometryGaiaReadStatus status =
            Internal::parseAstrometryGaiaPositions(input, gaia_xy, error);
        if (status == Internal::AstrometryGaiaReadStatus::Malformed) {
            MPIFailure::abortWorld(
                "parse PSF Gaia astrometry matches", filename + " " + error);
        }
        return gaia_xy;
    }

    // ==========================================
    // Structure: Store the shared minChi survivors for every exposure chip
    // Method: Address candidates by their original per-chip indices.
    // ==========================================
    using ActiveIndicesByChip = std::vector<std::vector<int>>;

    // ==========================================
    // Function: Build exposure-thresholded same-chip minChi survivor lists
    // Method: Select capped exposure-wide large-size references, compute every
    //         same-chip locus pair once, and threshold from reference-all pairs.
    // ==========================================
    static ActiveIndicesByChip buildMinChiActiveIndices(
        int nchip,
        ExposurePSFState& state) {
        std::vector<Internal::MinChiReferenceCandidate> reference_candidates;
        int locus_count = 0;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            const ChipPSFState& chip = state.chips[chip_index];
            for (int star_index = 0;
                 star_index < state.getNStar(chip_index); ++star_index) {
                if (!chip.selection[star_index].in_size_locus) continue;
                locus_count++;
                reference_candidates.push_back({
                    chip_index, star_index,
                    state.getStarPara(chip_index, star_index, 7)});
            }
        }

        const std::vector<Internal::MinChiReferenceCandidate> references =
            Internal::selectMinChiReferenceStars(
                reference_candidates,
                LensingConfig::psf_minchi_reference_fraction,
                LensingConfig::psf_minchi_reference_max_per_chip);
        std::vector<std::vector<bool>> is_reference(static_cast<std::size_t>(nchip));
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            is_reference[chip_index].assign(
                static_cast<std::size_t>(state.getNStar(chip_index)), false);
        }
        for (const Internal::MinChiReferenceCandidate& reference : references) {
            if (reference.chip_index >= 0 && reference.chip_index < nchip
                && reference.star_index >= 0
                && reference.star_index < state.getNStar(reference.chip_index)) {
                is_reference[reference.chip_index][reference.star_index] = true;
            }
        }

        std::vector<float> threshold_pair_chi;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            std::vector<Internal::MinChiCandidateView> candidates;
            candidates.reserve(chip.selection.size());
            for (int star_index = 0;
                 star_index < state.getNStar(chip_index); ++star_index) {
                candidates.push_back({
                    &chip.selection[star_index].chi_window,
                    chip.selection[star_index].in_size_locus,
                    is_reference[chip_index][star_index]});
            }
            Internal::MinChiPairResult pair_result =
                Internal::computeMinChiAndThresholdPairs(candidates);
            for (int star_index = 0;
                 star_index < state.getNStar(chip_index); ++star_index) {
                chip.selection[star_index].min_chi = pair_result.min_chi[star_index];
            }
            threshold_pair_chi.insert(
                threshold_pair_chi.end(), pair_result.threshold_pair_chi.begin(),
                pair_result.threshold_pair_chi.end());
        }

        const float min_chi_threshold = estimateUpperTailThreshold(
            threshold_pair_chi, LensingConfig::psf_minchi_sigma_cut);
        std::cout << "PSF_MINCHI exposure locus=" << locus_count
                  << " references=" << references.size()
                  << " threshold_pairs=" << threshold_pair_chi.size()
                  << " threshold=" << min_chi_threshold << std::endl;
        ActiveIndicesByChip active_indices(static_cast<std::size_t>(nchip));
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            const ChipPSFState& chip = state.chips[chip_index];
            for (int star_index = 0;
                 star_index < state.getNStar(chip_index); ++star_index) {
                const Internal::StarSelectionState& selection =
                    chip.selection[star_index];
                if (selection.in_size_locus && std::isfinite(selection.min_chi)
                    && selection.min_chi <= min_chi_threshold) {
                    active_indices[chip_index].push_back(star_index);
                }
            }
        }
        return active_indices;
    }

    // ==========================================
    // Function: Commit one direct Type-3 pre-PRESS survivor collection
    // Method: Validate original indices, enforce the chip minimum atomically,
    //         update legacy flags, and release every rejected temporary cache.
    // ==========================================
    static void commitAdaptivePairSelection(
        int nchip,
        const ActiveIndicesByChip& proposed_indices,
        ExposurePSFState& state) {
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            const int nstar = state.getNStar(chip_index);
            std::vector<bool> proposed(static_cast<std::size_t>(nstar), false);
            if (chip_index < static_cast<int>(proposed_indices.size())) {
                for (int star_index : proposed_indices[chip_index]) {
                    if (star_index >= 0 && star_index < nstar) {
                        proposed[star_index] = true;
                    }
                }
            }
            const int proposed_count = static_cast<int>(std::count(
                proposed.begin(), proposed.end(), true));
            const bool keep_chip =
                proposed_count >= LensingConfig::nstar_min_local;
            for (int star_index = 0; star_index < nstar; ++star_index) {
                Internal::StarSelectionState& selection =
                    chip.selection[star_index];
                selection.selected_group = keep_chip && proposed[star_index];
                state.getStarPara(chip_index, star_index, 4) =
                    selection.selected_group ? 1.0 : -1.0;
                if (!selection.selected_group) {
                    std::vector<float>().swap(selection.chi_window);
                }
            }
            std::cout << "PSF_TYPE3_CHIP chip=" << (chip_index + 1)
                      << " proposed=" << proposed_count
                      << " retained=" << (keep_chip ? proposed_count : 0)
                      << " decision="
                      << (keep_chip ? "KEEP" : "REJECT_MINIMUM")
                      << std::endl;
        }
    }

    // ==========================================
    // Function: Log one Type-3 adaptive histogram decision
    // Method: Publish FD sample/grid/topology fields for reproducible fail-open
    //         diagnosis without requiring a new rendered diagnostic product.
    // ==========================================
    static void logAdaptiveHistogram(
        const char* label,
        const Internal::PSFUpperElbowHistogramResult& result,
        const char* decision) {
        std::cout << label
                  << " finite=" << result.finite_value_count
                  << " fd_samples=" << result.fd_sample_count
                  << " fd_scale_samples="
                  << result.fd_scale_sample_count
                  << " iqr=" << result.fd_iqr
                  << " width=" << result.bin_width
                  << " origin=" << result.bin_origin
                  << " bins=" << result.histogram.size()
                  << " main_peak=" << result.main_peak_bin
                  << " rightmost_valid="
                  << result.rightmost_valid_peak_bin
                  << " first_invalid=" << result.first_invalid_peak_bin
                  << " elbow=" << result.elbow_bin
                  << " cut=" << result.cut
                  << " status="
                  << Internal::psfUpperElbowStatusName(result.status)
                  << " decision=" << decision << std::endl;
    }

    // ==========================================
    // Function: Select minChi survivors with adaptive pair/fraction elbows
    // Method: Estimate an exposure pair-chi cut, recompute same-chip finite
    //         pairs into per-star bad fractions, then commit a direct pre-PRESS gate.
    // ==========================================
    static void applyAdaptivePairFractionSelection(
        int nchip,
        ExposurePSFState& state,
        const ActiveIndicesByChip& active_indices) {
        for (ChipPSFState& chip : state.chips) {
            for (Internal::StarSelectionState& selection : chip.selection) {
                selection.bad_pair_fraction = 0.0;
            }
        }

        // ==========================================
        // Function: Derive the Type-3 Stage-1 FD scale count
        // Method: Sum minChi survivors with checked size_t addition; a
        //         representational overflow preserves the fail-open selection.
        // ==========================================
        std::size_t active_star_count = 0U;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            const std::size_t chip_active_count =
                active_indices[chip_index].size();
            if (chip_active_count
                > std::numeric_limits<std::size_t>::max()
                    - active_star_count) {
                Internal::PSFUpperElbowHistogramResult count_failure;
                count_failure.status =
                    Internal::PSFUpperElbowStatus::InvalidInput;
                logAdaptiveHistogram(
                    "PSF_TYPE3_PAIR", count_failure, "FAIL_OPEN");
                commitAdaptivePairSelection(nchip, active_indices, state);
                return;
            }
            active_star_count += chip_active_count;
        }

        std::vector<float> pair_chi;
        try {
            for (int chip_index = 0; chip_index < nchip; ++chip_index) {
                const ChipPSFState& chip = state.chips[chip_index];
                const std::vector<int>& active = active_indices[chip_index];
                for (std::size_t first = 0;
                     first + 1 < active.size(); ++first) {
                    for (std::size_t second = first + 1;
                         second < active.size(); ++second) {
                        const float chi = Internal::normalizedChiDistance(
                            chip.selection[active[first]].chi_window,
                            chip.selection[active[second]].chi_window);
                        if (std::isfinite(chi)) pair_chi.push_back(chi);
                    }
                }
            }
        } catch (const std::bad_alloc&) {
            Internal::PSFUpperElbowHistogramResult allocation_failure;
            allocation_failure.status =
                Internal::PSFUpperElbowStatus::AllocationFailure;
            logAdaptiveHistogram(
                "PSF_TYPE3_PAIR", allocation_failure, "FAIL_OPEN");
            commitAdaptivePairSelection(nchip, active_indices, state);
            return;
        } catch (const std::length_error&) {
            Internal::PSFUpperElbowHistogramResult allocation_failure;
            allocation_failure.status =
                Internal::PSFUpperElbowStatus::AllocationFailure;
            logAdaptiveHistogram(
                "PSF_TYPE3_PAIR", allocation_failure, "FAIL_OPEN");
            commitAdaptivePairSelection(nchip, active_indices, state);
            return;
        }

        Internal::PSFUpperElbowHistogramResult pair_result;
        const Internal::PSFUpperElbowHistogramConfig pair_config = {
            LensingConfig::psf_pair_chi_valid_peak_fraction,
            false,
            false,
            false,
            active_star_count,
            LensingConfig::psf_type3_elbow_search_height_fraction};
        if (!Internal::estimatePSFUpperElbowCut(
                pair_chi, pair_config, pair_result)) {
            logAdaptiveHistogram(
                "PSF_TYPE3_PAIR", pair_result, "FAIL_OPEN");
            commitAdaptivePairSelection(nchip, active_indices, state);
            return;
        }
        logAdaptiveHistogram("PSF_TYPE3_PAIR", pair_result, "APPLY");
        std::vector<float>().swap(pair_chi);

        std::vector<std::vector<std::size_t>> total_pairs;
        std::vector<std::vector<std::size_t>> bad_pairs;
        try {
            total_pairs.resize(static_cast<std::size_t>(nchip));
            bad_pairs.resize(static_cast<std::size_t>(nchip));
            for (int chip_index = 0; chip_index < nchip; ++chip_index) {
                total_pairs[chip_index].assign(
                    static_cast<std::size_t>(state.getNStar(chip_index)), 0U);
                bad_pairs[chip_index].assign(
                    static_cast<std::size_t>(state.getNStar(chip_index)), 0U);
                const ChipPSFState& chip = state.chips[chip_index];
                const std::vector<int>& active = active_indices[chip_index];
                for (std::size_t first = 0;
                     first + 1 < active.size(); ++first) {
                    for (std::size_t second = first + 1;
                         second < active.size(); ++second) {
                        const int first_index = active[first];
                        const int second_index = active[second];
                        const float chi = Internal::normalizedChiDistance(
                            chip.selection[first_index].chi_window,
                            chip.selection[second_index].chi_window);
                        if (!std::isfinite(chi)) continue;
                        total_pairs[chip_index][first_index]++;
                        total_pairs[chip_index][second_index]++;
                        if (Internal::isPSFType3BadPair(
                                static_cast<double>(chi), pair_result.cut)) {
                            bad_pairs[chip_index][first_index]++;
                            bad_pairs[chip_index][second_index]++;
                        }
                    }
                }
            }
        } catch (const std::bad_alloc&) {
            Internal::PSFUpperElbowHistogramResult allocation_failure;
            allocation_failure.status =
                Internal::PSFUpperElbowStatus::AllocationFailure;
            logAdaptiveHistogram(
                "PSF_TYPE3_FRACTION", allocation_failure, "FAIL_OPEN");
            commitAdaptivePairSelection(nchip, active_indices, state);
            return;
        } catch (const std::length_error&) {
            Internal::PSFUpperElbowHistogramResult allocation_failure;
            allocation_failure.status =
                Internal::PSFUpperElbowStatus::AllocationFailure;
            logAdaptiveHistogram(
                "PSF_TYPE3_FRACTION", allocation_failure, "FAIL_OPEN");
            commitAdaptivePairSelection(nchip, active_indices, state);
            return;
        }

        std::vector<double> fraction_values;
        std::size_t positive_fraction_count = 0;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            for (int star_index : active_indices[chip_index]) {
                const std::size_t denominator =
                    total_pairs[chip_index][star_index];
                if (denominator == 0U) continue;
                const double fraction = static_cast<double>(
                    bad_pairs[chip_index][star_index])
                    / static_cast<double>(denominator);
                chip.selection[star_index].bad_pair_fraction = fraction;
                fraction_values.push_back(fraction);
                if (fraction > 0.0) positive_fraction_count++;
            }
        }

        bool apply_fraction_cut = false;
        Internal::PSFUpperElbowHistogramResult fraction_result;
        if (positive_fraction_count == 0U) {
            fraction_result.finite_value_count = fraction_values.size();
            fraction_result.status =
                Internal::PSFUpperElbowStatus::NoFDSamples;
            logAdaptiveHistogram(
                "PSF_TYPE3_FRACTION",
                fraction_result,
                fraction_values.empty()
                    ? "NO_DENOMINATORS"
                    : "ALL_ZERO_PASS");
        } else {
            const Internal::PSFUpperElbowHistogramConfig fraction_config = {
                LensingConfig::psf_bad_fraction_valid_peak_fraction,
                true,
                true,
                true,
                0U,
                LensingConfig::psf_type3_elbow_search_height_fraction};
            apply_fraction_cut = Internal::estimatePSFUpperElbowCut(
                fraction_values, fraction_config, fraction_result);
            logAdaptiveHistogram(
                "PSF_TYPE3_FRACTION",
                fraction_result,
                apply_fraction_cut ? "APPLY" : "FAIL_OPEN");
        }

        ActiveIndicesByChip proposed_indices(static_cast<std::size_t>(nchip));
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            const std::vector<int>& active = active_indices[chip_index];
            std::vector<double> fractions(active.size(), 0.0);
            std::vector<bool> has_denominator(active.size(), false);
            for (std::size_t index = 0; index < active.size(); ++index) {
                const int star_index = active[index];
                has_denominator[index] =
                    total_pairs[chip_index][star_index] > 0U;
                fractions[index] =
                    state.chips[chip_index]
                        .selection[star_index].bad_pair_fraction;
            }
            const Internal::PSFType3ChipSelection chip_selection =
                Internal::selectPSFType3FractionSurvivors(
                    fractions,
                    has_denominator,
                    apply_fraction_cut,
                    fraction_result.cut,
                    LensingConfig::nstar_min_local);
            double fraction_sum = 0.0;
            double fraction_maximum = 0.0;
            for (std::size_t index = 0; index < active.size(); ++index) {
                if (!has_denominator[index]) continue;
                fraction_sum += fractions[index];
                fraction_maximum = std::max(
                    fraction_maximum, fractions[index]);
            }
            const double fraction_mean = chip_selection.finite_pair_count > 0U
                ? fraction_sum
                    / static_cast<double>(chip_selection.finite_pair_count)
                : 0.0;
            for (std::size_t index = 0;
                 index < chip_selection.selected.size(); ++index) {
                if (chip_selection.selected[index]) {
                    proposed_indices[chip_index].push_back(active[index]);
                }
            }
            std::cout << "PSF_TYPE3_FRACTION_CHIP chip="
                      << (chip_index + 1)
                      << " active=" << active.size()
                      << " finite_denominator="
                      << chip_selection.finite_pair_count
                      << " mean=" << fraction_mean
                      << " max=" << fraction_maximum
                      << " selected=" << chip_selection.retained_count
                      << " minimum_reject="
                      << (chip_selection.rejected_by_minimum ? 1 : 0)
                      << " decision="
                      << (chip_selection.rejected_by_minimum
                              ? "REJECT_MINIMUM"
                              : (chip_selection.finite_pair_count == 0U
                                     ? "NO_FINITE_PAIRS"
                                     : "KEEP"))
                      << std::endl;
        }
        commitAdaptivePairSelection(nchip, proposed_indices, state);
    }

    // ==========================================
    // Function: Select PSF stars from quality-valid candidates
    // Method: Apply common Gaia/star-area/minChi selection, dispatch grouping and the
    //         shared policy, then publish the pre-PRESS selected distribution.
    // ==========================================
    void starSelection(
        int nchip,
        const std::vector<std::string>& imageFiles,
        const std::string& dirOutput,
        ExposurePSFState& state) {
        int quality_valid_count = 0;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            for (int star_index = 0; star_index < state.getNStar(chip_index); ++star_index) {
                if (state.getStarPara(chip_index, star_index, 4) > 0.0) {
                    quality_valid_count++;
                }
            }
        }
        if (quality_valid_count < LensingConfig::psf_exposure_min_candidates) {
            rejectExposureCandidates(state);
            return;
        }

        std::vector<Internal::PSFCountSample> count_samples;
        count_samples.reserve(static_cast<std::size_t>(quality_valid_count));
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            const std::vector<std::array<double, 2>> gaia_xy =
                readAstrometryGaiaPositions(imageFiles[chip_index], dirOutput);
            for (int star_index = 0; star_index < state.getNStar(chip_index); ++star_index) {
                if (state.getStarPara(chip_index, star_index, 4) <= 0.0) continue;
                Internal::StarSelectionState& selection = chip.selection[star_index];
                selection.gaia_matched = Internal::hasNearestGaiaMatch(
                    state.getStarPara(chip_index, star_index, 1),
                    state.getStarPara(chip_index, star_index, 2),
                    gaia_xy,
                    LensingConfig::psf_gaia_match_radius_pix);
                const int star_area = static_cast<int>(std::llround(
                    state.getStarPara(
                        chip_index,
                        star_index,
                        ChipPSFState::star_area_index)));
                count_samples.push_back({
                    star_area,
                    selection.gaia_matched});
            }
        }

        Internal::PSFCountLocus locus;
        Internal::PSFCountLocusDiagnostics locus_diagnostics;
        const Internal::PSFCountLocusConfig locus_config = {
            LensingConfig::psf_count_pilot_clip_sigma,
            LensingConfig::psf_count_pilot_clip_iterations,
            LensingConfig::psf_count_zero_mad_quantile,
            LensingConfig::psf_count_hist_range_sigma,
            LensingConfig::psf_count_locus_sigma,
            LensingConfig::psf_count_locus_min_samples,
            LensingConfig::psf_gaia_locus_min_matches};
        if (!Internal::estimatePSFCountLocus(
                count_samples,
                locus_config,
                locus,
                &locus_diagnostics)) {
            rejectExposureCandidates(state);
            return;
        }

        const std::string prefix_e =
            UniversalUtils::getPrefixExpo(imageFiles[0]);
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            for (int star_index = 0; star_index < state.getNStar(chip_index); ++star_index) {
                Internal::StarSelectionState& selection = chip.selection[star_index];
                selection.selected_group = false;
                selection.selected_press = false;
                selection.min_chi = std::numeric_limits<float>::infinity();
                selection.bad_pair_fraction = 0.0;
                if (state.getStarPara(chip_index, star_index, 4) <= 0.0) continue;

                const double star_area = state.getStarPara(
                    chip_index,
                    star_index,
                    ChipPSFState::star_area_index);
                selection.in_size_locus =
                    star_area > locus.lower && star_area < locus.upper;
                if (!selection.in_size_locus
                    || selection.full_power_sum <= 0.0
                    || selection.chi_window.empty()) {
                    state.getStarPara(chip_index, star_index, 4) = -1.0;
                    std::vector<float>().swap(selection.chi_window);
                    continue;
                }
                const double inverse_sum = 1.0 / selection.full_power_sum;
                for (float& value : selection.chi_window) {
                    value = static_cast<float>(static_cast<double>(value) * inverse_sum);
                }
            }
        }

        const ActiveIndicesByChip active_indices =
            buildMinChiActiveIndices(nchip, state);
        std::size_t minchi_survivor_total = 0;
        for (const std::vector<int>& chip_indices : active_indices) {
            minchi_survivor_total += chip_indices.size();
        }
        std::vector<int> minchi_survivor_star_areas;
        minchi_survivor_star_areas.reserve(minchi_survivor_total);
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            for (int star_index : active_indices[chip_index]) {
                minchi_survivor_star_areas.push_back(
                    static_cast<int>(std::llround(
                        state.getStarPara(
                            chip_index,
                            star_index,
                            ChipPSFState::star_area_index))));
            }
        }
        Internal::populateMinChiSurvivorCountHistogram(
            minchi_survivor_star_areas, locus_diagnostics);
        applyAdaptivePairFractionSelection(nchip, state, active_indices);

        std::vector<int> selected_group_star_areas;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            const ChipPSFState& chip = state.chips[chip_index];
            for (int star_index = 0;
                 star_index < state.getNStar(chip_index); ++star_index) {
                if (chip.selection[star_index].selected_group) {
                    selected_group_star_areas.push_back(
                        static_cast<int>(std::llround(
                            state.getStarPara(
                                chip_index,
                                star_index,
                                ChipPSFState::star_area_index))));
                }
            }
        }
        Internal::populateSelectedGroupCountHistogram(
            selected_group_star_areas, locus_diagnostics);
        writePSFCountLocusSVG(
            dirOutput, prefix_e, locus, locus_diagnostics);
    }

    // ==========================================
    // Structure: Hold one chip's ordered polynomial fitting samples
    // Method: Preserve original candidate indices alongside positions and
    //         contiguous row-major power stamps.
    // ==========================================
    struct ChipFitSamples {
        std::vector<int> star_indices;
        std::vector<std::array<double, 2>> positions;
        std::vector<float> power;
    };

    // ==========================================
    // Function: Read all candidate power stamps for one chip
    // Method: Validate and return the direct catalog-shaped Stage-4 stamp cube.
    // ==========================================
    static std::vector<float> readChipCandidatePower(
        int chip_index,
        const std::vector<std::string>& imageFiles,
        const std::string& dirOutput,
        const ExposurePSFState& state) {
        const int nstar = state.getNStar(chip_index);
        if (nstar <= 0) return {};
        const int ns = LensingConfig::ns;
        const std::string prefix =
            UniversalUtils::getPrefix(imageFiles[chip_index]);
        const std::string filename = OutputLayout::chipPath(
            dirOutput, "stamps/fits_StarCanP", prefix,
            "_star_can_power.fits");
        std::vector<float> power;
        readExpectedStampCube(
            filename, ns, ns, nstar, power,
            "read PRESS PSF star-power cube");
        return power;
    }

    // ==========================================
    // Function: Gather an ordered subset of one chip's fitting samples
    // Method: Validate each requested original candidate and copy its position
    //         and full power stamp into a compact fitting buffer.
    // ==========================================
    static bool buildChipFitSamples(
        int chip_index,
        const std::vector<int>& selected_indices,
        const std::vector<float>& all_power,
        const ExposurePSFState& state,
        ChipFitSamples& samples) {
        samples = {};
        const int ns = LensingConfig::ns;
        const int pixel_count = ns * ns;
        samples.star_indices.reserve(selected_indices.size());
        samples.positions.reserve(selected_indices.size());
        samples.power.reserve(selected_indices.size()
                              * static_cast<std::size_t>(pixel_count));
        for (int star_index : selected_indices) {
            if (star_index < 0 || star_index >= state.getNStar(chip_index)) return false;
            const double x = state.getStarPara(chip_index, star_index, 1);
            const double y = state.getStarPara(chip_index, star_index, 2);
            const double size = state.getStarPara(chip_index, star_index, 7);
            const double e1 = state.getStarPara(chip_index, star_index, 8);
            const double e2 = state.getStarPara(chip_index, star_index, 9);
            const std::size_t offset =
                static_cast<std::size_t>(star_index) * pixel_count;
            if (!isFinitePSFStar(
                    all_power, offset, pixel_count, x, y, size, e1, e2)) {
                return false;
            }
            samples.star_indices.push_back(star_index);
            samples.positions.push_back({x, y});
            samples.power.insert(
                samples.power.end(),
                all_power.begin() + static_cast<std::ptrdiff_t>(offset),
                all_power.begin() + static_cast<std::ptrdiff_t>(offset + pixel_count));
        }
        return true;
    }

    // ==========================================
    // Function: Fit one ordered chip sample set and compute leverage
    // Method: Use Lite's sole normalized local-polynomial design with runtime
    //         chip geometry and share coefficients/hat diagonals with outputs.
    // ==========================================
    static LinearSolve::SolveStatus fitChipSamples(
        const ChipFitSamples& samples,
        std::vector<double>& coefficients,
        std::vector<double>& leverage,
        LinearSolve::SolveDiagnostics& diagnostics) {
        const int sample_count = static_cast<int>(samples.star_indices.size());
        const LensingRuntimeConfig& lensing =
            RuntimeConfigStore::get().lensing;
        return itpNormPSF(
            sample_count, samples.power, samples.positions,
            LensingConfig::ns, LensingConfig::npl,
            lensing.chipnx, lensing.chipny,
            coefficients, &diagnostics, &leverage);
    }

    // ==========================================
    // Function: Mark one chip invalid after a PRESS fitting failure
    // Method: Clear cached fit products and reject every group-selected candidate.
    // ==========================================
    static void invalidatePressChip(int chip_index, ExposurePSFState& state) {
        ChipPSFState& chip = state.chips[chip_index];
        chip.fit.clear();
        for (int star_index = 0; star_index < state.getNStar(chip_index); ++star_index) {
            if (chip.selection[star_index].selected_group) {
                state.getStarPara(chip_index, star_index, 4) = -1.0;
                std::vector<float>().swap(
                    chip.selection[star_index].chi_window);
            }
            chip.selection[star_index].selected_press = false;
        }
    }

    // ==========================================
    // Function: Name one non-mutating PRESS removal decision
    // Method: Map the pure safeguard result to a stable diagnostic label.
    // ==========================================
    static const char* pressRemovalDecisionName(
        Internal::PressRemovalDecision decision) {
        switch (decision) {
            case Internal::PressRemovalDecision::Disabled: return "DISABLED";
            case Internal::PressRemovalDecision::NoOutliers: return "NO_OUTLIER";
            case Internal::PressRemovalDecision::TooManyOutliers: return "TOO_MANY_OUTLIERS";
            case Internal::PressRemovalDecision::WouldUnderrunMinimum: return "WOULD_UNDERRUN_MINIMUM";
            case Internal::PressRemovalDecision::Apply: return "APPLY";
        }
        return "UNKNOWN";
    }

    // ==========================================
    // Function: Apply optional exposure-wide standardized PRESS cleanup
    // Method: Preserve a valid first fit, propose removals from leverage-corrected
    //         scores, guard them, refit in temporary state, and commit on success.
    // ==========================================
    void applyPressSelection(
        int nchip,
        const std::vector<std::string>& imageFiles,
        const std::string& dirOutput,
        ExposurePSFState& state) {
        const LensingRuntimeConfig& lensing =
            RuntimeConfigStore::get().lensing;
        const int ns = LensingConfig::ns;
        const int pixel_count = ns * ns;
        const Internal::PSFChiWindow chi_window =
            Internal::getPSFChiWindow(ns);
        std::vector<float> exposure_standardized_scores;

        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            chip.fit.clear();
            std::vector<int> selected_indices;
            for (int star_index = 0; star_index < state.getNStar(chip_index); ++star_index) {
                chip.selection[star_index].selected_press = false;
                chip.selection[star_index].press_raw_score = 0.0;
                chip.selection[star_index].press_standardized_score = 0.0;
                chip.selection[star_index].leverage = 0.0;
                if (chip.selection[star_index].selected_group) {
                    selected_indices.push_back(star_index);
                }
            }
            if (static_cast<int>(selected_indices.size())
                < LensingConfig::nstar_min_local) {
                invalidatePressChip(chip_index, state);
                continue;
            }

            const std::vector<float> all_power = readChipCandidatePower(
                chip_index, imageFiles, dirOutput, state);
            ChipFitSamples samples;
            if (!buildChipFitSamples(
                    chip_index, selected_indices, all_power, state, samples)) {
                LinearSolve::reportFailure(
                    "PSFModel::applyPressSelection",
                    LinearSolve::SolveStatus::FailedSolver,
                    "chip=" + std::to_string(chip_index + 1)
                        + " reason=INVALID_SELECTED_SAMPLE action=MARK_CHIP_INVALID");
                invalidatePressChip(chip_index, state);
                continue;
            }

            LinearSolve::SolveDiagnostics diagnostics;
            std::vector<double> coefficients;
            std::vector<double> leverage;
            const LinearSolve::SolveStatus status = fitChipSamples(
                samples, coefficients, leverage, diagnostics);
            if (status != LinearSolve::SolveStatus::Normal) {
                LinearSolve::reportFailure(
                    "PSFModel::applyPressSelection", status,
                    "chip=" + std::to_string(chip_index + 1) + " "
                        + LinearSolve::diagnosticsContext(diagnostics)
                        + " action=MARK_CHIP_INVALID");
                invalidatePressChip(chip_index, state);
                continue;
            }

            bool loo_valid = true;
            std::vector<float> chip_standardized_scores;
            chip_standardized_scores.reserve(samples.star_indices.size());
            for (int local_index = 0;
                 local_index < static_cast<int>(samples.star_indices.size());
                 ++local_index) {
                const double x = 2.0 * (
                    samples.positions[local_index][0]
                    / static_cast<double>(lensing.chipnx)) - 1.0;
                const double y = 2.0 * (
                    samples.positions[local_index][1]
                    / static_cast<double>(lensing.chipny)) - 1.0;
                std::vector<float> model;
                std::vector<float> constant_model;
                getPSFModel(
                    ns, LensingConfig::npl, coefficients,
                    x, y, model, constant_model);
                double squared_sum = 0.0;
                for (int row = chi_window.first;
                     row <= chi_window.last && loo_valid; ++row) {
                    for (int column = chi_window.first;
                         column <= chi_window.last; ++column) {
                        const int pixel = row * ns + column;
                        const double observed = samples.power[
                            static_cast<std::size_t>(local_index) * pixel_count + pixel];
                        double loo_residual = 0.0;
                        double loo_model = 0.0;
                        if (!Internal::computeAnalyticLOO(
                                observed, model[pixel], leverage[local_index],
                                LensingConfig::psf_loo_min_denom,
                                loo_residual, loo_model)) {
                            loo_valid = false;
                            break;
                        }
                        squared_sum += loo_residual * loo_residual;
                    }
                }
                if (!loo_valid) break;
                const double raw_press_score = std::sqrt(
                    squared_sum / static_cast<double>(chi_window.pixelCount()));
                double standardized_press_score = 0.0;
                if (!Internal::computeLeverageStandardizedPress(
                        raw_press_score, leverage[local_index],
                        LensingConfig::psf_loo_min_denom,
                        standardized_press_score)) {
                    loo_valid = false;
                    break;
                }
                const int original_index = samples.star_indices[local_index];
                chip.selection[original_index].press_raw_score = raw_press_score;
                chip.selection[original_index].press_standardized_score = standardized_press_score;
                chip.selection[original_index].leverage = leverage[local_index];
                chip_standardized_scores.push_back(static_cast<float>(standardized_press_score));
            }
            if (!loo_valid) {
                LinearSolve::reportFailure(
                    "PSFModel::applyPressSelection",
                    LinearSolve::SolveStatus::FailedIllConditioned,
                    "chip=" + std::to_string(chip_index + 1)
                        + " reason=INVALID_LOO action=MARK_CHIP_INVALID");
                invalidatePressChip(chip_index, state);
                continue;
            }

            exposure_standardized_scores.insert(
                exposure_standardized_scores.end(),
                chip_standardized_scores.begin(), chip_standardized_scores.end());

            chip.fit.valid = true;
            chip.fit.initial_star_count = static_cast<int>(samples.star_indices.size());
            chip.fit.star_indices = samples.star_indices;
            chip.fit.coefficients = std::move(coefficients);
            chip.fit.leverage = std::move(leverage);
            for (int star_index : chip.fit.star_indices) {
                chip.selection[star_index].selected_press = true;
                state.getStarPara(chip_index, star_index, 4) = 1.0;
            }
        }

        const float press_threshold = estimateUpperTailThreshold(
            exposure_standardized_scores, LensingConfig::psf_press_sigma_cut);
        std::cout << "PSF_PRESS exposure standardized_scores="
                  << exposure_standardized_scores.size()
                  << " threshold=" << press_threshold << std::endl;
        for (int chip_index = 0; chip_index < nchip; ++chip_index) {
            ChipPSFState& chip = state.chips[chip_index];
            if (!chip.fit.valid) continue;

            const std::vector<int> first_fit_indices = chip.fit.star_indices;
            std::vector<int> rejected_indices;
            std::vector<int> retained_indices;
            retained_indices.reserve(chip.fit.star_indices.size());
            for (int star_index : chip.fit.star_indices) {
                const bool retained = std::isfinite(
                        chip.selection[star_index].press_standardized_score)
                    && chip.selection[star_index].press_standardized_score <= press_threshold;
                if (retained) {
                    retained_indices.push_back(star_index);
                } else {
                    rejected_indices.push_back(star_index);
                }
            }

            const Internal::PressRemovalDecision decision = Internal::decidePressRemoval(
                LensingConfig::psf_press_rejection_enabled,
                static_cast<int>(first_fit_indices.size()),
                static_cast<int>(rejected_indices.size()),
                LensingConfig::nstar_min_local,
                LensingConfig::psf_press_max_removals);
            if (decision != Internal::PressRemovalDecision::Apply) {
                chip.fit.press_removed_any = false;
                std::cout << "PSF_PRESS chip=" << (chip_index + 1)
                          << " first_fit=" << first_fit_indices.size()
                          << " flagged=" << rejected_indices.size()
                          << " decision=" << pressRemovalDecisionName(decision)
                          << " final=" << first_fit_indices.size() << std::endl;
                continue;
            }

            const std::vector<float> all_power = readChipCandidatePower(
                chip_index, imageFiles, dirOutput, state);
            ChipFitSamples retained_samples;
            if (!buildChipFitSamples(
                    chip_index, retained_indices, all_power, state,
                    retained_samples)) {
                LinearSolve::reportFailure(
                    "PSFModel::applyPressSelectionRefit",
                    LinearSolve::SolveStatus::FailedSolver,
                    "chip=" + std::to_string(chip_index + 1)
                        + " reason=INVALID_RETAINED_SAMPLE action=PRESS_REFIT_FALLBACK");
                chip.fit.press_removed_any = false;
                std::cout << "PSF_PRESS chip=" << (chip_index + 1)
                          << " first_fit=" << first_fit_indices.size()
                          << " flagged=" << rejected_indices.size()
                          << " decision=PRESS_REFIT_FALLBACK final="
                          << first_fit_indices.size() << std::endl;
                continue;
            }

            LinearSolve::SolveDiagnostics diagnostics;
            std::vector<double> coefficients;
            std::vector<double> leverage;
            const LinearSolve::SolveStatus status = fitChipSamples(
                retained_samples, coefficients, leverage, diagnostics);
            if (status != LinearSolve::SolveStatus::Normal) {
                LinearSolve::reportFailure(
                    "PSFModel::applyPressSelectionRefit", status,
                    "chip=" + std::to_string(chip_index + 1) + " "
                        + LinearSolve::diagnosticsContext(diagnostics)
                        + " action=PRESS_REFIT_FALLBACK");
                chip.fit.press_removed_any = false;
                std::cout << "PSF_PRESS chip=" << (chip_index + 1)
                          << " first_fit=" << first_fit_indices.size()
                          << " flagged=" << rejected_indices.size()
                          << " decision=PRESS_REFIT_FALLBACK final="
                          << first_fit_indices.size() << std::endl;
                continue;
            }

            if (!chip.fit.tryCommitPressRefit(
                    true, retained_samples.star_indices,
                    std::move(coefficients), std::move(leverage))) {
                LinearSolve::reportFailure(
                    "PSFModel::applyPressSelectionRefit",
                    LinearSolve::SolveStatus::FailedSolver,
                    "chip=" + std::to_string(chip_index + 1)
                        + " reason=INVALID_REFIT_CACHE action=PRESS_REFIT_FALLBACK");
                chip.fit.press_removed_any = false;
                continue;
            }
            std::vector<bool> retained_mask(
                static_cast<std::size_t>(state.getNStar(chip_index)), false);
            for (int star_index : chip.fit.star_indices) retained_mask[star_index] = true;
            for (int star_index : first_fit_indices) {
                const bool retained = retained_mask[star_index];
                chip.selection[star_index].selected_press = retained;
                state.getStarPara(chip_index, star_index, 4) = retained ? 1.0 : -1.0;
                if (!retained) std::vector<float>().swap(chip.selection[star_index].chi_window);
            }
            for (std::size_t local_index = 0;
                 local_index < chip.fit.star_indices.size(); ++local_index) {
                const int original_index = chip.fit.star_indices[local_index];
                chip.selection[original_index].leverage = chip.fit.leverage[local_index];
            }
            std::cout << "PSF_PRESS chip=" << (chip_index + 1)
                      << " first_fit=" << first_fit_indices.size()
                      << " flagged=" << rejected_indices.size()
                      << " decision=REMOVAL_APPLIED final="
                      << chip.fit.star_indices.size() << std::endl;
        }
    }

    // ==========================================
    // Function: Assemble exposure-wide selected-star Fourier power
    // Method: Append each actual-size chip cube, retain the original selection
    //         order, and publish a dense exposure-level cube without padding.
    // ==========================================
    void plotStarExpo(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        int ns = LensingConfig::ns;

        std::vector<int> opt;
        std::vector<float> star_test;

        int ntot = 0;
        int start = 0;

        for (int ichip = 0; ichip < nchip; ++ichip) {
            const int nstar = state.getNStar(ichip);
            if (nstar == 0) continue;

            std::string prefix = UniversalUtils::getPrefix(imageFiles[ichip]);
            std::string filepath = OutputLayout::chipPath(
                dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");

            std::vector<float> star;
            readExpectedStampCube(
                filepath, ns, ns, nstar, star,
                "read exposure PSF star-power cube");

            star_test.insert(star_test.end(), star.begin(), star.end());
            opt.resize(static_cast<std::size_t>(start + nstar), 0);

            for (int i = 0; i < nstar; ++i) {
                if (state.getStarPara(ichip, i, 4) <= 0.0) continue;
                ntot++;
                opt[start + i] = 1;
            }
            start += nstar;
        }

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string out_filename = dirOutput + "/stamps/fits_StarP/" + prefix_e + "_star_power_expo.fits";

        if (ntot > 0) {
            FitsIO::writeSelectedStampCube(
                out_filename, ns, ns, star_test, start, opt, 1);
        }
    }

    // ==========================================
    // Function: Publish exposure-level star diagnostics and PSF visualization
    // Method: Aggregate selected stars and route text/FITS products through
    //         checked main-process writers.
    // ==========================================
    void plotStars(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, int nc, const std::vector<std::array<double, 4>>& p_chip, ExposurePSFState& state) {
        (void)nc;
        int nm = 1000;
        int nstar_min_local = LensingConfig::nstar_min_local;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        // prefix_dir inlined: per-type stamps/ subdirs (reorganized layout)

        std::string info_filename = dirOutput + "/stamps/dat_StarInfo/" + prefix_e + "_star_info_expo.dat";
        MainIO::OutputFile outfile(info_filename);
        outfile << std::setprecision(10);

        outfile << "# ichip nstar FWHM e1 e2 chi_d\n";

        const std::size_t initial_star_capacity = checkedElementCount(
            {static_cast<std::size_t>(
                 RuntimeConfigStore::get().lensing.nmax_chip),
             static_cast<std::size_t>(LensingConfig::nstar_max)},
            "reserve exposure PSF diagnostics");
        std::vector<std::array<double, 5>> sk;
        sk.reserve(initial_star_capacity);

        for (int k = 0; k < nchip; ++k) {
            double FWHM_ave = 0.0;
            double e1_ave = 0.0;
            double e2_ave = 0.0;
            double chi_d_ave = 0.0;
            int nums = 0;
            int prev_idx = -1;

            for (int i = 0; i < state.getNStar(k); ++i) {
                if (state.getStarPara(k, i, 4) <= 0.0) continue;
                nums++;
                sk.push_back({
                    state.getStarPara(k, i, 5),
                    state.getStarPara(k, i, 6),
                    state.getStarPara(k, i, 7),
                    state.getStarPara(k, i, 8),
                    state.getStarPara(k, i, 9)
                });

                FWHM_ave += state.getStarPara(k, i, 10);
                e1_ave += state.getStarPara(k, i, 8);
                e2_ave += state.getStarPara(k, i, 9);
                if (nums >= 2 && prev_idx != -1) {
                    chi_d_ave += Internal::normalizedChiDistance(
                        state.chips[k].selection[i].chi_window,
                        state.chips[k].selection[prev_idx].chi_window);
                }
                prev_idx = i;
            }

            if (nums >= nstar_min_local) {
                FWHM_ave /= nums;
                e1_ave /= nums;
                e2_ave /= nums;
                chi_d_ave /= (nums - 1.0);
                outfile << (k + 1) << " " << nums << " "
                        << std::scientific << std::setprecision(10)
                        << FWHM_ave << " " << e1_ave << " " << e2_ave << " " << chi_d_ave << "\n";
            } else {
                outfile << (k + 1) << " 0 -99.0 -99.0 -99.0 -99.0\n";
            }
        }

        // std::cout << prefix_dir << " total no. of stars: " << sk.size() << std::endl;
        outfile.close();

        std::vector<float> PSFmap;
        ImageProcessing::drawShearExpo(nm, PSFmap, p_chip, sk, 200.0, 1.0);

        std::string fits_filename = dirOutput + "/stamps/fits_PsfSrc/" + prefix_e + "_PSF_source.fits";
        FitsIO::writeImage(fits_filename, nm, nm, PSFmap);
    }

    // ==========================================
    // Function: Fit and serialize local PSF models.
    // Method: Preserve F77 model layout with 17-digit double serialization.
    // ==========================================
    void makePSFLocalFit(int nchip, const std::vector<std::string>& imageFiles, const std::string& dirOutput, ExposurePSFState& state) {
        int ns = LensingConfig::ns;
        int npl = LensingConfig::npl;
        const LensingRuntimeConfig& lensing =
            RuntimeConfigStore::get().lensing;

        std::string prefix_e = UniversalUtils::getPrefixExpo(imageFiles[0]);
        std::string comp_filename = dirOutput + "/stamps/dat_StarComp/" + prefix_e + "_star_comp_expo.dat";
        MainIO::OutputFile file90(comp_filename);
        if (!file90.is_open()) {
            std::cerr << "makePSFLocalFit: Error opening " << comp_filename << std::endl;
            return;
        }
        file90 << std::setprecision(17);

        for (int k = 0; k < nchip; ++k) {
            int nums = 0;
            std::string prefix = UniversalUtils::getPrefix(imageFiles[k]);

            std::vector<float> star;
            if (state.getNStar(k) > 0) {
                std::string filepath = OutputLayout::chipPath(
                    dirOutput, "stamps/fits_StarCanP", prefix, "_star_can_power.fits");
                readExpectedStampCube(
                    filepath, ns, ns, state.getNStar(k), star,
                    "read local-fit PSF star-power cube");
            }

            std::string coe_filename = OutputLayout::chipPath(
                dirOutput, "stamps/dat_PsfFit", prefix, "_PSF_coe_local.dat");
            MainIO::OutputFile file10(coe_filename);
            if (!file10.is_open()) {
                std::cerr << "makePSFLocalFit: Error opening " << coe_filename << std::endl;
                continue;
            }
            file10 << std::setprecision(17);

            std::vector<std::array<double, 2>> posi;
            std::vector<std::array<double, 3>> sshape;
            std::vector<float> star_local;
            int removed_non_finite = 0;

            const Internal::ChipPSFFitState& cached_fit = state.chips[k].fit;
            for (int i : cached_fit.star_indices) {
                double px = state.getStarPara(k, i, 1);
                double py = state.getStarPara(k, i, 2);
                double shape_size = state.getStarPara(k, i, 7);
                double shape_e1 = state.getStarPara(k, i, 8);
                double shape_e2 = state.getStarPara(k, i, 9);
                std::size_t star_offset = static_cast<std::size_t>(i) * ns * ns;
                if (!isFinitePSFStar(star, star_offset, ns * ns, px, py,
                                     shape_size, shape_e1, shape_e2)) {
                    removed_non_finite++;
                    continue;
                }
                nums++;
                posi.push_back({px, py});
                sshape.push_back({shape_size, shape_e1, shape_e2});
                for (int idx = 0; idx < ns * ns; ++idx) {
                    star_local.push_back(star[star_offset + idx]);
                }
            }
            std::vector<double> PSF_coe_l = cached_fit.coefficients;
            std::vector<double> final_leverage = cached_fit.leverage;
            LinearSolve::SolveDiagnostics fit_diagnostics;
            LinearSolve::SolveStatus fit_status = cached_fit.valid
                ? LinearSolve::SolveStatus::Normal
                : LinearSolve::SolveStatus::FailedRankDeficient;

            if (nums >= LensingConfig::nstar_min_local &&
                fit_status == LinearSolve::SolveStatus::Normal &&
                final_leverage.size() == static_cast<std::size_t>(nums)) {

                file90 << (k + 1) << " " << nums << " 1\n";

                std::vector<float> poly_cochi2(nums);
                float poly_ave = 0.0f, poly_std = 0.0f;

                for (int i = 0; i < nums; ++i) {
                    double xx = 2.0 * (posi[i][0] / static_cast<double>(lensing.chipnx)) - 1.0;
                    double yy = 2.0 * (posi[i][1] / static_cast<double>(lensing.chipny)) - 1.0;
                    std::vector<float> model, model0;
                    getPSFModel(ns, npl, PSF_coe_l, xx, yy, model, model0);
                    ExStar::anaChi2Simple(ns, model.data(), model0.data(), poly_cochi2[i]);

                    std::vector<float> loo_model(static_cast<std::size_t>(ns) * ns);
                    for (int idx = 0; idx < ns * ns; ++idx) {
                        const double observed = star_local[
                            static_cast<std::size_t>(i) * ns * ns + idx];
                        double residual_value = 0.0;
                        double model_value = 0.0;
                        if (!Internal::computeAnalyticLOO(
                                observed, model[idx], final_leverage[i],
                                LensingConfig::psf_loo_min_denom,
                                residual_value, model_value)) {
                            MPIFailure::abortWorld(
                                "generate final Lite PSF LOO diagnostics",
                                "exposure=" + prefix_e
                                    + " chip=" + std::to_string(k + 1)
                                    + " star=" + std::to_string(i));
                        }
                        loo_model[idx] = static_cast<float>(model_value);
                    }

                    std::array<double, 2> ee = {0.0, 0.0};
                    double size = 0.0;
                    getPowerAll(ns, ns, loo_model, ee, size, 0.02f);

                    double msshape_size = size;
                    double msshape_e1 = ee[0];
                    double msshape_e2 = ee[1];

                    float px = static_cast<float>(posi[i][0]);
                    float py = static_cast<float>(posi[i][1]);

                    file90 << px << " " << py << " "
                           << sshape[i][0] << " " << sshape[i][1] << " " << sshape[i][2] << " "
                           << msshape_size << " " << msshape_e1 << " " << msshape_e2 << "\n";
                }

                ExStar::getArrayAveStd(poly_cochi2, poly_ave, poly_std);
                file10 << nums << " 1 " << poly_ave << " " << poly_std << "\n";
                for (int i = 0; i < ns; ++i) {
                    for (int j = 0; j < ns; ++j) {
                        for (int u = 0; u < npl + 1; ++u) {
                            file10 << PSF_coe_l[(j * ns + i) * (npl + 1) + u] << (u == npl ? "" : " ");
                        }
                        file10 << "\n";
                    }
                }
            } else {
                if (nums < LensingConfig::nstar_min_local) {
                    LinearSolve::reportFailure(
                        "PSFModel::itpNormPSF", LinearSolve::SolveStatus::FailedRankDeficient,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " valid_samples=" + std::to_string(nums) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " required=" + std::to_string(LensingConfig::nstar_min_local) +
                            " action=MARK_CHIP_INVALID");
                } else {
                    LinearSolve::reportFailure(
                        "PSFModel::itpNormPSF", fit_status,
                        "exposure=" + prefix_e + " chip=" + std::to_string(k + 1) +
                            " " + LinearSolve::diagnosticsContext(fit_diagnostics) +
                            " removed_samples=" + std::to_string(removed_non_finite) +
                            " action=MARK_CHIP_INVALID");
                }
                file10 << nums << " -1 -1 -1\n";
                file90 << (k + 1) << " " << nums << " -1\n";
            }
            file10.close();
        }
        file90.close();
    }

    // Mathematical Interpolation Helpers

    // ==========================================
    // Function: Fit all PSF-frequency pixels with a shared spatial design
    // Method: Factor the constant and polynomial designs once, then reuse both QR objects for every RHS pixel.
    // Note:   Lite retains only the normalized local-polynomial design, so the
    //         chip-coordinate mapping to [-1,1] is unconditional.
    // ==========================================
    static LinearSolve::SolveStatus fitPSFCoefficients(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny,
        std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics,
        std::vector<double>* leverage) {
        const int pixel_count = ns * ns;
        PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);

        LinearSolve::SolveDiagnostics local_diagnostics;
        LinearSolve::SolveDiagnostics& diag = diagnostics == nullptr ? local_diagnostics : *diagnostics;
        diag = {};
        diag.rows = nsam;
        diag.cols = npp;
        diag.required_rank = npp;

        if (nsam < npp || ns <= 0 || npp <= 0 ||
            static_cast<int>(posi.size()) < nsam ||
            image.size() < static_cast<size_t>(nsam) * pixel_count ||
            nx <= 0 || ny <= 0) {
            diag.rank = std::min(nsam, npp);
            return LinearSolve::SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd polynomial_design(nsam, npp);
        Eigen::MatrixXd constant_design = Eigen::MatrixXd::Ones(nsam, 1);
        for (int sample = 0; sample < nsam; ++sample) {
            double xx = 2.0 * (posi[sample][0] / static_cast<double>(nx)) - 1.0;
            double yy = 2.0 * (posi[sample][1] / static_cast<double>(ny)) - 1.0;
            if (!std::isfinite(xx) || !std::isfinite(yy)) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
            for (int term = 0; term < npp; ++term) {
                polynomial_design(sample, term) = UniversalUtils::fitFunc2(xx, yy, term);
            }
        }
        if (!polynomial_design.allFinite()) {
            return LinearSolve::SolveStatus::FailedSolver;
        }
        for (std::size_t idx = 0; idx < static_cast<std::size_t>(nsam) * pixel_count; ++idx) {
            if (!std::isfinite(image[idx])) {
                return LinearSolve::SolveStatus::FailedSolver;
            }
        }

        LinearSolve::LeastSquaresQR constant_solver;
        LinearSolve::SolveDiagnostics constant_diagnostics;
        LinearSolve::SolveStatus status = constant_solver.factorize(constant_design, constant_diagnostics);
        if (status != LinearSolve::SolveStatus::Normal) {
            diag = constant_diagnostics;
            return status;
        }

        LinearSolve::LeastSquaresQR polynomial_solver;
        status = polynomial_solver.factorize(polynomial_design, diag);
        if (status != LinearSolve::SolveStatus::Normal) {
            return status;
        }

        if (leverage != nullptr) {
            Eigen::MatrixXd covariance_base;
            status = polynomial_solver.unscaledCovariance(covariance_base);
            if (status != LinearSolve::SolveStatus::Normal) return status;
            leverage->assign(static_cast<std::size_t>(nsam), 0.0);
            for (int sample = 0; sample < nsam; ++sample) {
                const double value = (
                    polynomial_design.row(sample) * covariance_base
                    * polynomial_design.row(sample).transpose())(0, 0);
                double loo_residual = 0.0;
                double loo_model = 0.0;
                if (!Internal::computeAnalyticLOO(
                        0.0, 0.0, value,
                        LensingConfig::psf_loo_min_denom,
                        loo_residual, loo_model)) {
                    leverage->clear();
                    return LinearSolve::SolveStatus::FailedIllConditioned;
                }
                (*leverage)[sample] = std::max(0.0, value);
            }
        }

        Eigen::VectorXd rhs(nsam);
        Eigen::VectorXd constant_solution;
        Eigen::VectorXd polynomial_solution;
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                const int pixel = j * ns + i;
                for (int sample = 0; sample < nsam; ++sample) {
                    rhs(sample) = image[static_cast<size_t>(sample) * pixel_count + pixel];
                }

                status = constant_solver.solve(rhs, constant_solution);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);
                    return status;
                }
                status = polynomial_solver.solve(rhs, polynomial_solution);
                if (status != LinearSolve::SolveStatus::Normal) {
                    PSF_coe.assign(static_cast<size_t>(pixel_count) * (npp + 1), 0.0);
                    return status;
                }

                for (int term = 0; term < npp; ++term) {
                    PSF_coe[static_cast<size_t>(pixel) * (npp + 1) + term] = polynomial_solution(term);
                }
                PSF_coe[static_cast<size_t>(pixel) * (npp + 1) + npp] = constant_solution(0);
            }
        }
        return LinearSolve::SolveStatus::Normal;
    }

    // ==========================================
    // Function: Fit normalized local PSF coefficients
    // Method: Normalize chip coordinates to [-1,1] and delegate to the shared-design batch fitter.
    // ==========================================
    LinearSolve::SolveStatus itpNormPSF(
        int nsam, const std::vector<float>& image,
        const std::vector<std::array<double, 2>>& posi,
        int ns, int npp, int nx, int ny, std::vector<double>& PSF_coe,
        LinearSolve::SolveDiagnostics* diagnostics,
        std::vector<double>* leverage) {
        return fitPSFCoefficients(
            nsam, image, posi, ns, npp, nx, ny,
            PSF_coe, diagnostics, leverage);
    }

    void getPSFModel(int ns, int npp, const std::vector<double>& PSF_coe, double xx, double yy, std::vector<float>& modelp, std::vector<float>& model0) {
        modelp.assign(static_cast<size_t>(ns) * ns, 0.0f);
        model0.assign(static_cast<size_t>(ns) * ns, 0.0f);
        std::vector<double> coep(npp);
        std::vector<double> coe0(1);
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                for (int k = 0; k < npp; ++k) {
                    coep[k] = PSF_coe[(j * ns + i) * (npp + 1) + k];
                }
                coe0[0] = PSF_coe[(j * ns + i) * (npp + 1) + npp];
                model0[j * ns + i] = static_cast<float>(UniversalUtils::funcVal2(xx, yy, 1, coe0));
                modelp[j * ns + i] = static_cast<float>(UniversalUtils::funcVal2(xx, yy, npp, coep));
                if (std::isnan(modelp[j * ns + i])) {
                    modelp[0] = modelp[j * ns + i];
                    return;
                }
            }
        }
    }

    void getPowerArea(int nx, int ny, const std::vector<float>& power, int& area, float thresh_ratio) {
        int cx = nx / 2;
        int cy = ny / 2;
        float thresh = power[cy * nx + cx] * thresh_ratio;

        std::vector<int> mark(nx * ny, 0);
        std::vector<int> stack_x(nx * ny, 0);
        std::vector<int> stack_y(nx * ny, 0);

        int area_cnt = 0;
        int area0 = 0;

        mark[cy * nx + cx] = 1;
        stack_x[0] = cx;
        stack_y[0] = cy;
        area_cnt = 1;

        while (area_cnt > area0) {
            int tempi = area_cnt;
            for (int i = area0; i < tempi; ++i) {
                int x = stack_x[i];
                int y = stack_y[i];
                for (int u = std::max(x - 1, 0); u <= std::min(x + 1, nx - 1); ++u) {
                    for (int v = std::max(y - 1, 0); v <= std::min(y + 1, ny - 1); ++v) {
                        int idx = v * nx + u;
                        if (mark[idx] == 0 && power[idx] >= thresh) {
                            mark[idx] = 1;
                            stack_x[area_cnt] = u;
                            stack_y[area_cnt] = v;
                            area_cnt++;
                        }
                    }
                }
            }
            area0 = tempi;
        }

        area = (area_cnt - 1) / 2;
    }

    void getPowerE(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, float thresh_ratio) {
        int cx = nx / 2;
        int cy = ny / 2;
        float thresh = power[cy * nx + cx] * thresh_ratio;

        e[0] = 0.0;
        e[1] = 0.0;
        double norm = 0.0;

        std::vector<int> mark(nx * ny, 0);
        std::vector<int> stack_x(nx * ny, 0);
        std::vector<int> stack_y(nx * ny, 0);

        int area_cnt = 0;
        int area0 = 0;

        mark[cy * nx + cx] = 1;
        stack_x[0] = cx;
        stack_y[0] = cy;
        area_cnt = 1;

        while (area_cnt > area0) {
            int tempi = area_cnt;
            for (int i = area0; i < tempi; ++i) {
                int x = stack_x[i];
                int y = stack_y[i];
                for (int u = std::max(x - 1, 0); u <= std::min(x + 1, nx - 1); ++u) {
                    for (int v = std::max(y - 1, 0); v <= std::min(y + 1, ny - 1); ++v) {
                        int idx = v * nx + u;
                        if (mark[idx] == 0 && power[idx] >= thresh) {
                            mark[idx] = 1;
                            stack_x[area_cnt] = u;
                            stack_y[area_cnt] = v;
                            area_cnt++;

                            double kx = u - cx;
                            double ky = v - cy;
                            double p_val = power[idx];
                            e[0] += p_val * (kx * kx - ky * ky);
                            e[1] += p_val * 2.0 * kx * ky;
                            norm += p_val * (kx * kx + ky * ky);
                        }
                    }
                }
            }
            area0 = tempi;
        }

        e[0] /= norm;
        e[1] /= norm;
    }

    void getPowerAll(int nx, int ny, const std::vector<float>& power, std::array<double, 2>& e, double& size, float thresh_ratio) {
        int area = 0;
        getPowerArea(nx, ny, power, area, thresh_ratio);
        size = area;
        getPowerE(nx, ny, power, e, thresh_ratio);
    }

    // ==========================================
    // Function: Measure historical PSF FWHM and exact exp(-1) star area
    // Method: Count integer threshold pixels once, then reuse the unchanged
    //         area-minus-1e-5 conversion for the reported FWHM.
    // ==========================================
    void getPSFFWHM(
        const std::vector<float>& power,
        double& FWHM,
        int& star_area) {
        star_area = Internal::countPSFStarArea(power, LensingConfig::ns);
        FWHM = Internal::fwhmFromStarArea(
            static_cast<double>(star_area),
            LensingConfig::ns,
            LensingConfig::pixel_size);
    }
}
