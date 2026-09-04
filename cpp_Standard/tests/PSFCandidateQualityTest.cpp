#include "process_main/PSFCandidateQuality.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using PSFModel::Internal::CandidatePowerStatus;
using PSFModel::Internal::assessCandidatePower;
using PSFModel::Internal::candidateDiagnosticsAreFinite;

// ==========================================
// Function: Stop the candidate-quality test on one failed invariant
// Method: Print a focused message and return a nonzero process status.
// ==========================================
void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "PSF candidate quality test failed: " << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

// ==========================================
// Function: Build a small accepted corrected-power spectrum
// Method: Give DC and its eight neighbours a positive core while retaining
//         one allowed negative outer Fourier pixel.
// ==========================================
std::vector<float> makeAcceptedPower() {
    std::vector<float> power(25, 0.0f);
    const int neighbors[] = {6, 7, 8, 11, 13, 16, 17, 18};
    power[12] = 10.0f;
    for (int index : neighbors) power[index] = 1.0f;
    power[0] = -1.0f;
    return power;
}

// ==========================================
// Function: Verify accepted and boundary spectra
// Method: Check signed-sum reuse and the inclusive zero-median boundary.
// ==========================================
void testAcceptedPower() {
    double sum_power = 0.0;
    double chi_window_sum = 0.0;
    std::vector<float> power = makeAcceptedPower();
    require(assessCandidatePower(5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::Accepted,
            "finite positive spectrum must pass");
    require(std::abs(sum_power - 17.0) < 1.0e-12,
            "accepted spectrum must return its signed sum");
    require(chi_window_sum > 0.0,
            "accepted spectrum must return a positive chi-window sum");

    const int neighbors[] = {6, 7, 8, 11, 13, 16, 17, 18};
    for (int i = 0; i < 4; ++i) power[neighbors[i]] = -1.0f;
    for (int i = 4; i < 8; ++i) power[neighbors[i]] = 1.0f;
    power[0] = 0.0f;
    require(assessCandidatePower(5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::Accepted,
            "zero central-neighbour median must remain accepted");
}

// ==========================================
// Function: Verify non-finite corrected-power rejection
// Method: Inject NaN and both infinity signs into an otherwise valid spectrum.
// ==========================================
void testNonFinitePower() {
    const float bad_values[] = {
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity(),
        -std::numeric_limits<float>::infinity()
    };
    for (float bad_value : bad_values) {
        std::vector<float> power = makeAcceptedPower();
        power[3] = bad_value;
        double sum_power = 0.0;
        double chi_window_sum = 0.0;
        require(assessCandidatePower(
                    5, 5, power, sum_power, chi_window_sum)
                    == CandidatePowerStatus::NonFinitePower,
                "every NaN/Inf sign must be rejected");
    }
}

// ==========================================
// Function: Verify negative low-frequency-core rejection
// Method: Make five of eight DC neighbours negative while keeping total power positive.
// ==========================================
void testNegativeCoreMedian() {
    std::vector<float> power(25, 0.0f);
    const int neighbors[] = {6, 7, 8, 11, 13, 16, 17, 18};
    power[12] = 20.0f;
    for (int i = 0; i < 5; ++i) power[neighbors[i]] = -1.0f;
    for (int i = 5; i < 8; ++i) power[neighbors[i]] = 1.0f;
    double sum_power = 0.0;
    double chi_window_sum = 0.0;
    require(assessCandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::NegativeCoreMedian,
            "negative eight-neighbour median must be rejected");
}

// ==========================================
// Function: Verify non-positive normalization rejection
// Method: Preserve a positive central core while forcing signed total power to zero and negative.
// ==========================================
void testNonPositiveSum() {
    std::vector<float> power = makeAcceptedPower();
    double sum_power = 0.0;
    double chi_window_sum = 0.0;
    power[0] = -18.0f;
    require(assessCandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::NonPositiveSum,
            "zero signed sum must be rejected");
    power[0] = -19.0f;
    require(assessCandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::NonPositiveSum,
            "negative signed sum must be rejected");
}

// ==========================================
// Function: Verify non-positive central chi-window rejection
// Method: Keep the full signed sum positive while forcing the shared 5x5-test
//         window sum first to zero and then negative.
// ==========================================
void testNonPositiveChiWindowSum() {
    std::vector<float> power = makeAcceptedPower();
    double sum_power = 0.0;
    double chi_window_sum = 0.0;
    power[0] = -13.0f;
    require(assessCandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::NonPositiveChiWindowSum,
            "zero chi-window sum must be rejected after a positive full sum");
    require(sum_power > 0.0 && chi_window_sum == 0.0,
            "zero-window fixture must preserve the intended sums");

    power[0] = -14.0f;
    require(assessCandidatePower(
                5, 5, power, sum_power, chi_window_sum)
                == CandidatePowerStatus::NonPositiveChiWindowSum,
            "negative chi-window sum must be rejected after a positive full sum");
    require(sum_power > 0.0 && chi_window_sum < 0.0,
            "negative-window fixture must preserve the intended sums");
}

// ==========================================
// Function: Verify structural and derived-diagnostic guards
// Method: Reject incomplete images and every non-finite diagnostic field.
// ==========================================
void testStructuralAndDiagnosticValidity() {
    double sum_power = 0.0;
    double chi_window_sum = 0.0;
    require(assessCandidatePower(
                5, 5, std::vector<float>(24), sum_power, chi_window_sum)
                == CandidatePowerStatus::InvalidShape,
            "incomplete power image must be rejected");
    require(assessCandidatePower(
                2, 2, std::vector<float>(4), sum_power, chi_window_sum)
                == CandidatePowerStatus::InvalidShape,
            "image without a complete DC neighbourhood must be rejected");
    require(candidateDiagnosticsAreFinite(1.0, 0.1, -0.1),
            "finite diagnostics must pass");
    require(!candidateDiagnosticsAreFinite(
                1.0, std::numeric_limits<double>::quiet_NaN(), 0.0),
            "NaN ellipticity must be rejected before Stage 8");
    require(!candidateDiagnosticsAreFinite(
                std::numeric_limits<double>::infinity(), 0.0, 0.0),
            "infinite size must be rejected before Stage 8");
}

}  // namespace

// ==========================================
// Function: Run the corrected-power candidate quality regression suite
// Method: Exercise all gate outcomes and report one success line.
// ==========================================
int main() {
    testAcceptedPower();
    testNonFinitePower();
    testNegativeCoreMedian();
    testNonPositiveSum();
    testNonPositiveChiWindowSum();
    testStructuralAndDiagnosticValidity();
    std::cout << "PSF candidate quality tests passed\n";
    return EXIT_SUCCESS;
}
