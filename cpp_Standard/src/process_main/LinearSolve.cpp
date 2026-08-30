#include "process_main/LinearSolve.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace LinearSolve {

    namespace {

        // ==========================================
        // Function: Compute a dense matrix one-norm
        // Method: Return the maximum absolute column sum for a small QR triangular factor.
        // ==========================================
        double matrixOneNorm(const Eigen::MatrixXd& matrix) {
            return matrix.cwiseAbs().colwise().sum().maxCoeff();
        }

        // ==========================================
        // Function: Compute a dense matrix infinity-norm
        // Method: Return the maximum absolute row sum for a small QR triangular factor.
        // ==========================================
        double matrixInfinityNorm(const Eigen::MatrixXd& matrix) {
            return matrix.cwiseAbs().rowwise().sum().maxCoeff();
        }

    }

    // ==========================================
    // Function: Factor and classify a least-squares design matrix
    // Method: Apply column scaling and pivoted QR, prove clearly healthy cases with a conservative R-factor bound,
    //         otherwise compute exact singular values of R; reject every non-NORMAL condition class.
    // ==========================================
    SolveStatus LeastSquaresQR::factorize(const Eigen::MatrixXd& design, SolveDiagnostics& diagnostics) {
        rows_ = design.rows();
        cols_ = design.cols();
        ready_ = false;

        diagnostics = {};
        diagnostics.rows = rows_;
        diagnostics.cols = cols_;
        diagnostics.required_rank = cols_;
        diagnostics.warning_rcond_threshold = condition_warning_rcond;
        diagnostics.rank_tolerance =
            static_cast<double>(std::max(rows_, cols_)) *
            std::numeric_limits<double>::epsilon();
        diagnostics.severe_rcond_threshold = std::min(
            diagnostics.warning_rcond_threshold,
            std::max(condition_severe_rcond_floor,
                     condition_severe_rank_multiplier * diagnostics.rank_tolerance));

        if (!design.allFinite()) {
            return SolveStatus::FailedSolver;
        }
        if (rows_ <= 0 || cols_ <= 0 || rows_ < cols_) {
            diagnostics.rank = std::min(rows_, cols_);
            diagnostics.condition_number = std::numeric_limits<double>::infinity();
            diagnostics.condition_class = ConditionClass::RankDeficient;
            return SolveStatus::FailedRankDeficient;
        }

        column_scales_.resize(cols_);
        Eigen::MatrixXd scaled_design = design;
        for (Eigen::Index column = 0; column < cols_; ++column) {
            double scale = design.col(column).stableNorm();
            if (!std::isfinite(scale)) {
                return SolveStatus::FailedSolver;
            }
            if (scale == 0.0) {
                diagnostics.rank = 0;
                diagnostics.condition_number = std::numeric_limits<double>::infinity();
                diagnostics.condition_class = ConditionClass::RankDeficient;
                return SolveStatus::FailedRankDeficient;
            }
            column_scales_(column) = scale;
            scaled_design.col(column) /= scale;
        }

        qr_.compute(scaled_design);
        qr_.setThreshold(diagnostics.rank_tolerance);

        if (!qr_.matrixQR().allFinite()) {
            return SolveStatus::FailedSolver;
        }

        diagnostics.rank = qr_.rank();
        if (diagnostics.rank < diagnostics.required_rank) {
            diagnostics.condition_number = std::numeric_limits<double>::infinity();
            diagnostics.condition_class = ConditionClass::RankDeficient;
            return SolveStatus::FailedRankDeficient;
        }

        Eigen::MatrixXd upper = qr_.matrixR().topLeftCorner(cols_, cols_)
                                    .template triangularView<Eigen::Upper>();
        Eigen::MatrixXd inverse_upper =
            upper.template triangularView<Eigen::Upper>().solve(
                Eigen::MatrixXd::Identity(cols_, cols_));
        if (inverse_upper.allFinite()) {
            const double condition_one =
                matrixOneNorm(upper) * matrixOneNorm(inverse_upper);
            const double condition_infinity =
                matrixInfinityNorm(upper) * matrixInfinityNorm(inverse_upper);
            const double condition_upper_bound =
                std::sqrt(condition_one * condition_infinity);
            if (std::isfinite(condition_upper_bound) && condition_upper_bound > 0.0) {
                diagnostics.condition_number = condition_upper_bound;
                diagnostics.reciprocal_condition = 1.0 / condition_upper_bound;
                if (diagnostics.reciprocal_condition >=
                    diagnostics.warning_rcond_threshold) {
                    diagnostics.estimated_reliable_digits = std::max(
                        0.0,
                        -std::log10(std::numeric_limits<double>::epsilon()) -
                            std::log10(diagnostics.condition_number));
                    diagnostics.condition_class = ConditionClass::Healthy;
                    ready_ = true;
                    return SolveStatus::Normal;
                }
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> singular_solver(upper);
        const Eigen::VectorXd singular_values = singular_solver.singularValues();
        diagnostics.condition_exact = true;
        if (singular_values.size() != cols_ || !singular_values.allFinite()) {
            return SolveStatus::FailedSolver;
        }

        diagnostics.sigma_max = singular_values(0);
        diagnostics.sigma_min = singular_values(cols_ - 1);
        if (!(diagnostics.sigma_max > 0.0) || !std::isfinite(diagnostics.sigma_max)) {
            diagnostics.condition_number = std::numeric_limits<double>::infinity();
            diagnostics.condition_class = ConditionClass::RankDeficient;
            return SolveStatus::FailedRankDeficient;
        }

        diagnostics.reciprocal_condition =
            diagnostics.sigma_min / diagnostics.sigma_max;
        diagnostics.condition_number = diagnostics.reciprocal_condition > 0.0
                                           ? 1.0 / diagnostics.reciprocal_condition
                                           : std::numeric_limits<double>::infinity();

        const double absolute_rank_threshold =
            diagnostics.rank_tolerance * diagnostics.sigma_max;
        diagnostics.rank = 0;
        for (Eigen::Index index = 0; index < singular_values.size(); ++index) {
            if (singular_values(index) > absolute_rank_threshold) {
                diagnostics.rank++;
            }
        }
        if (diagnostics.rank < diagnostics.required_rank) {
            diagnostics.condition_class = ConditionClass::RankDeficient;
            return SolveStatus::FailedRankDeficient;
        }

        diagnostics.estimated_reliable_digits = std::max(
            0.0,
            -std::log10(std::numeric_limits<double>::epsilon()) -
                std::log10(diagnostics.condition_number));

        if (!std::isfinite(diagnostics.reciprocal_condition) ||
            !std::isfinite(diagnostics.condition_number)) {
            return SolveStatus::FailedSolver;
        }
        if (diagnostics.reciprocal_condition < diagnostics.severe_rcond_threshold) {
            diagnostics.condition_class = ConditionClass::Severe;
            return SolveStatus::FailedIllConditioned;
        }
        if (diagnostics.reciprocal_condition < diagnostics.warning_rcond_threshold) {
            diagnostics.condition_class = ConditionClass::Warning;
            return SolveStatus::ValidWithWarning;
        }

        diagnostics.condition_class = ConditionClass::Healthy;
        ready_ = true;
        return SolveStatus::Normal;
    }

    // ==========================================
    // Function: Solve one least-squares RHS
    // Method: Use the retained QR once and reject non-finite input or output without a fallback solver.
    // ==========================================
    SolveStatus LeastSquaresQR::solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution) const {
        if (!ready_ || rhs.rows() != rows_) {
            return SolveStatus::FailedSolver;
        }
        if (!rhs.allFinite()) {
            return SolveStatus::FailedSolver;
        }

        Eigen::VectorXd scaled_solution = qr_.solve(rhs);
        solution = scaled_solution.array() / column_scales_.array();
        if (solution.rows() != cols_ || !solution.allFinite()) {
            return SolveStatus::FailedSolver;
        }
        return SolveStatus::Normal;
    }

    // ==========================================
    // Function: Solve multiple least-squares RHS columns
    // Method: Apply the retained QR to a finite RHS matrix and reject any non-finite result.
    // ==========================================
    SolveStatus LeastSquaresQR::solve(const Eigen::MatrixXd& rhs, Eigen::MatrixXd& solution) const {
        if (!ready_ || rhs.rows() != rows_) {
            return SolveStatus::FailedSolver;
        }
        if (!rhs.allFinite()) {
            return SolveStatus::FailedSolver;
        }

        Eigen::MatrixXd scaled_solution = qr_.solve(rhs);
        solution = column_scales_.cwiseInverse().asDiagonal() * scaled_solution;
        if (solution.rows() != cols_ || solution.cols() != rhs.cols() || !solution.allFinite()) {
            return SolveStatus::FailedSolver;
        }
        return SolveStatus::Normal;
    }

    // ==========================================
    // Function: Construct the unscaled coefficient covariance factor
    // Method: Form P R^-1 R^-T P^T from the full-rank pivoted QR without explicitly inverting A^T A.
    // ==========================================
    SolveStatus LeastSquaresQR::unscaledCovariance(Eigen::MatrixXd& covariance) const {
        if (!ready_) {
            return SolveStatus::FailedSolver;
        }

        Eigen::MatrixXd upper = qr_.matrixR().topLeftCorner(cols_, cols_)
                                    .template triangularView<Eigen::Upper>();
        Eigen::MatrixXd inverse_upper = Eigen::MatrixXd::Zero(cols_, cols_);
        for (Eigen::Index column = 0; column < cols_; ++column) {
            for (Eigen::Index row = cols_; row-- > 0;) {
                double diagonal = upper(row, row);
                if (!std::isfinite(diagonal) || diagonal == 0.0) {
                    return SolveStatus::FailedSolver;
                }
                double value = row == column ? 1.0 : 0.0;
                for (Eigen::Index inner = row + 1; inner < cols_; ++inner) {
                    value -= upper(row, inner) * inverse_upper(inner, column);
                }
                inverse_upper(row, column) = value / diagonal;
            }
        }
        if (!inverse_upper.allFinite()) {
            return SolveStatus::FailedSolver;
        }

        Eigen::MatrixXd pivoted_covariance = inverse_upper * inverse_upper.transpose();
        Eigen::MatrixXd scaled_covariance =
            qr_.colsPermutation() * pivoted_covariance * qr_.colsPermutation().transpose();
        Eigen::VectorXd inverse_scales = column_scales_.cwiseInverse();
        covariance = inverse_scales.asDiagonal() * scaled_covariance * inverse_scales.asDiagonal();
        if (!covariance.allFinite()) {
            return SolveStatus::FailedSolver;
        }
        return SolveStatus::Normal;
    }

    // ==========================================
    // Function: Convert a numerical status to text
    // Method: Keep failure tokens stable for stdout parsing.
    // ==========================================
    const char* statusName(SolveStatus status) {
        switch (status) {
            case SolveStatus::Normal:
                return "NORMAL";
            case SolveStatus::ValidWithWarning:
                return "VALID_WITH_WARNING";
            case SolveStatus::FailedRankDeficient:
                return "FAILED_RANK_DEFICIENT";
            case SolveStatus::FailedIllConditioned:
                return "FAILED_ILL_CONDITIONED";
            case SolveStatus::FailedSolver:
                return "FAILED_SOLVER";
            case SolveStatus::FailedNegativeCovariance:
                return "FAILED_NEGATIVE_COVARIANCE";
        }
        return "UNKNOWN_FAILURE";
    }

    // ==========================================
    // Function: Convert a matrix-condition class to text
    // Method: Keep condition tokens stable for diagnostics and regression tests.
    // ==========================================
    const char* conditionClassName(ConditionClass condition_class) {
        switch (condition_class) {
            case ConditionClass::NotEvaluated:
                return "NOT_EVALUATED";
            case ConditionClass::Healthy:
                return "HEALTHY";
            case ConditionClass::Warning:
                return "WARNING";
            case ConditionClass::Severe:
                return "SEVERE";
            case ConditionClass::RankDeficient:
                return "RANK_DEFICIENT";
        }
        return "NOT_EVALUATED";
    }

    // ==========================================
    // Function: Format rank and condition diagnostics
    // Method: Emit stable key-value fields at 17-digit precision for one standardized error line.
    // ==========================================
    std::string diagnosticsContext(const SolveDiagnostics& diagnostics) {
        std::ostringstream context;
        context << std::scientific << std::setprecision(17)
                << "rows=" << diagnostics.rows
                << " cols=" << diagnostics.cols
                << " rank=" << diagnostics.rank
                << " required=" << diagnostics.required_rank
                << " condition_class=" << conditionClassName(diagnostics.condition_class)
                << " rcond=" << diagnostics.reciprocal_condition
                << " condition_number=" << diagnostics.condition_number
                << " sigma_min=" << diagnostics.sigma_min
                << " sigma_max=" << diagnostics.sigma_max
                << " rank_tolerance=" << diagnostics.rank_tolerance
                << " severe_threshold=" << diagnostics.severe_rcond_threshold
                << " warning_threshold=" << diagnostics.warning_rcond_threshold;
        context << " condition_exact=" << (diagnostics.condition_exact ? 1 : 0);
        return context.str();
    }

    // ==========================================
    // Function: Classify a covariance eigenvalue spectrum
    // Method: Reject only values below the absolute negative threshold and count positive modes at machine precision.
    // ==========================================
    SolveStatus analyzeCovarianceSpectrum(
        const Eigen::VectorXd& eigenvalues, int sample_count,
        int requested_modes, double negative_threshold,
        EigenSpectrumDiagnostics& diagnostics) {
        diagnostics = {};
        if (sample_count < 2 || requested_modes < 0 || eigenvalues.size() <= 0) {
            return SolveStatus::FailedRankDeficient;
        }
        if (!eigenvalues.allFinite()) {
            return SolveStatus::FailedSolver;
        }

        diagnostics.lambda_min = eigenvalues(0);
        diagnostics.lambda_max = eigenvalues(eigenvalues.size() - 1);
        if (diagnostics.lambda_min < negative_threshold) {
            return SolveStatus::FailedNegativeCovariance;
        }

        if (diagnostics.lambda_max > 0.0) {
            diagnostics.positive_tolerance =
                static_cast<double>(std::max<Eigen::Index>(eigenvalues.size(), sample_count)) *
                std::numeric_limits<double>::epsilon() * diagnostics.lambda_max;
        }
        for (Eigen::Index index = 0; index < eigenvalues.size(); ++index) {
            if (eigenvalues(index) > diagnostics.positive_tolerance) {
                diagnostics.positive_modes++;
            }
        }
        diagnostics.effective_modes = std::min(
            requested_modes,
            std::min(sample_count - 1, diagnostics.positive_modes));
        return SolveStatus::Normal;
    }

    // ==========================================
    // Function: Print a standardized non-NORMAL numerical result
    // Method: Suppress NORMAL, then emit exactly one Error record for warnings and failures.
    // ==========================================
    void reportFailure(const std::string& location, SolveStatus status, const std::string& context) {
        if (status == SolveStatus::Normal) {
            return;
        }
        std::ostringstream message;
        message << "Error: (" << location << ", " << statusName(status) << ")";
        if (!context.empty()) {
            message << " " << context;
        }
        message << '\n';
        std::cout << message.str();
        std::cout.flush();
    }

}
