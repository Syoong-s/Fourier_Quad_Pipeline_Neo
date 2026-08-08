#ifndef LINEAR_SOLVE_HPP
#define LINEAR_SOLVE_HPP

#include <Eigen/Dense>
#include <cstddef>
#include <string>

namespace LinearSolve {

    enum class SolveStatus {
        Normal,
        ValidWithWarning,
        FailedRankDeficient,
        FailedIllConditioned,
        FailedSolver,
        FailedNegativeCovariance
    };

    enum class ConditionClass {
        NotEvaluated,
        Healthy,
        Warning,
        Severe,
        RankDeficient
    };

    inline constexpr double condition_warning_rcond = 1.0e-8;
    inline constexpr double condition_severe_rcond_floor = 1.0e-12;
    inline constexpr double condition_severe_rank_multiplier = 100.0;

    struct SolveDiagnostics {
        Eigen::Index rows = 0;
        Eigen::Index cols = 0;
        Eigen::Index rank = 0;
        Eigen::Index required_rank = 0;
        std::size_t removed_samples = 0;
        double sigma_max = 0.0;
        double sigma_min = 0.0;
        double reciprocal_condition = 0.0;
        double condition_number = 0.0;
        double rank_tolerance = 0.0;
        double severe_rcond_threshold = 0.0;
        double warning_rcond_threshold = condition_warning_rcond;
        double estimated_reliable_digits = 0.0;
        bool condition_exact = false;
        ConditionClass condition_class = ConditionClass::NotEvaluated;
    };

    struct EigenSpectrumDiagnostics {
        double lambda_min = 0.0;
        double lambda_max = 0.0;
        double positive_tolerance = 0.0;
        int positive_modes = 0;
        int effective_modes = 0;
    };

    // ==========================================
    // Class: Reusable rank-revealing least-squares factorization
    // Method: Factor a finite design matrix once with column-pivoted QR and reuse it for one or more RHS vectors.
    // ==========================================
    class LeastSquaresQR {
    public:
        SolveStatus factorize(const Eigen::MatrixXd& design, SolveDiagnostics& diagnostics);
        SolveStatus solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution) const;
        SolveStatus solve(const Eigen::MatrixXd& rhs, Eigen::MatrixXd& solution) const;
        SolveStatus unscaledCovariance(Eigen::MatrixXd& covariance) const;

    private:
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_;
        Eigen::VectorXd column_scales_;
        Eigen::Index rows_ = 0;
        Eigen::Index cols_ = 0;
        bool ready_ = false;
    };

    // ==========================================
    // Function: Convert a numerical status to the stable stdout token
    // Method: Return a fixed uppercase string for every status value.
    // ==========================================
    const char* statusName(SolveStatus status);

    // ==========================================
    // Function: Convert a matrix-condition class to a stable token
    // Method: Return one of NOT_EVALUATED, HEALTHY, WARNING, SEVERE, or RANK_DEFICIENT.
    // ==========================================
    const char* conditionClassName(ConditionClass condition_class);

    // ==========================================
    // Function: Format reusable numerical diagnostics
    // Method: Serialize rank and condition fields with full double precision for one stdout error context.
    // ==========================================
    std::string diagnosticsContext(const SolveDiagnostics& diagnostics);

    // ==========================================
    // Function: Classify a covariance eigenvalue spectrum
    // Method: Apply the approved absolute negative threshold and machine-positive rank threshold without rejecting expected zero modes.
    // ==========================================
    SolveStatus analyzeCovarianceSpectrum(
        const Eigen::VectorXd& eigenvalues, int sample_count,
        int requested_modes, double negative_threshold,
        EigenSpectrumDiagnostics& diagnostics);

    // ==========================================
    // Function: Emit one standardized non-NORMAL numerical record
    // Method: Suppress NORMAL and write warnings or failures as Error: (location, type) with optional context.
    // ==========================================
    void reportFailure(const std::string& location, SolveStatus status, const std::string& context = "");

}

#endif // LINEAR_SOLVE_HPP
