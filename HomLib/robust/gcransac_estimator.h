#pragma once
#ifdef HOMLIB_HAS_GCRANSAC

#include <cmath>
#include <vector>
#include <memory>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "estimators/estimator.h"
#include "model.h"
#include "posedata.hpp"
#include "gcransac_model_utils.h"
#include "solver_utils.h"

namespace HomLib {

// ---------------------------------------------------------------------------
// GCRANSAC estimator with distortion-aware residuals.
// Implements the gcransac::estimator::Estimator interface.
// Used by both GCRANSAC and (via HomLibMAGSACEstimator) MAGSAC.
// ---------------------------------------------------------------------------
template <class Solver,
          class _MinimalSolverEngine,
          class _NonMinimalSolverEngine>
class HomLibGCEstimator
    : public gcransac::estimator::Estimator<cv::Mat, gcransac::Model> {
protected:
    const std::shared_ptr<_MinimalSolverEngine> minimal_solver;
    const std::shared_ptr<_NonMinimalSolverEngine> non_minimal_solver;
    Solver solver_;

public:
    HomLibGCEstimator(const Solver& solver,
                      const NonMinimalSettings& nm = {},
                      const NonMinimalEstimateFn& nm_solver_fn = {})
        : minimal_solver(std::make_shared<_MinimalSolverEngine>(solver)),
          non_minimal_solver(std::make_shared<_NonMinimalSolverEngine>(solver, nm, nm_solver_fn)),
          solver_(solver) {}

    static constexpr size_t sampleSize() {
        return _MinimalSolverEngine::sampleSize();
    }

    static constexpr size_t nonMinimalSampleSize() {
        return _NonMinimalSolverEngine::sampleSize();
    }

    static constexpr bool isWeightingApplicable() { return true; }

    static constexpr size_t maximumMinimalSolutions() {
        return _MinimalSolverEngine::maximumSolutions();
    }

    OLGA_INLINE size_t inlierLimit() const override {
        return 7 * sampleSize();
    }

    const _MinimalSolverEngine* getMinimalSolver() const {
        return minimal_solver.get();
    }

    const _NonMinimalSolverEngine* getNonMinimalSolver() const {
        return non_minimal_solver.get();
    }

    _MinimalSolverEngine* getMutableMinimalSolver() {
        return minimal_solver.get();
    }

    _NonMinimalSolverEngine* getMutableNonMinimalSolver() {
        return non_minimal_solver.get();
    }

    // Minimal estimation
    OLGA_INLINE bool estimateModel(
        const cv::Mat& data_,
        const size_t* sample_,
        std::vector<gcransac::Model>* models_) const override
    {
        return minimal_solver->estimateModel(
            data_, sample_, sampleSize(), *models_);
    }

    // Non-minimal estimation (no external normalization — HomLib solvers
    // handle normalization internally)
    OLGA_INLINE bool estimateModelNonminimal(
        const cv::Mat& data_,
        const size_t* sample_,
        const size_t& sample_number_,
        std::vector<gcransac::Model>* models_,
        const double* weights_ = nullptr) const override
    {
        if (sample_number_ < nonMinimalSampleSize())
            return false;
        return non_minimal_solver->estimateModel(
            data_, sample_, sample_number_, *models_, weights_);
    }

    // Distortion-aware squared residual
    OLGA_INLINE double squaredResidual(
        const cv::Mat& point_,
        const gcransac::Model& model_) const override
    {
        const double* s = reinterpret_cast<const double*>(point_.data);
        PoseData pose = gc_model_to_pose(model_);
        Eigen::Vector2d x(s[0], s[1]);
        Eigen::Vector2d y(s[2], s[3]);
        return detail::compute_squared_residual(solver_, pose, x, y);
    }

    OLGA_INLINE double residual(
        const cv::Mat& point_,
        const gcransac::Model& model_) const override
    {
        return std::sqrt(squaredResidual(point_, model_));
    }

    OLGA_INLINE bool isValidModel(
        gcransac::Model& model_,
        const cv::Mat& data_,
        const std::vector<size_t>& inliers_,
        const size_t* minimal_sample_,
        const double threshold_,
        bool& model_updated_) const override
    {
        const double n = model_.descriptor.block<3, 3>(0, 0).norm();
        return n > 1e-15;
    }

    // Oriented constraint check (same as RobustHomographyEstimator)
    OLGA_INLINE bool isValidSample(
        const cv::Mat& data_,
        const size_t* sample_) const override
    {
        if constexpr (sampleSize() < 4)
            return true;

        const double* a = reinterpret_cast<const double*>(data_.row(sample_[0]).data);
        const double* b = reinterpret_cast<const double*>(data_.row(sample_[1]).data);
        const double* c = reinterpret_cast<const double*>(data_.row(sample_[2]).data);
        const double* d = reinterpret_cast<const double*>(data_.row(sample_[3]).data);

        Eigen::Vector3d p, q;

        p << a[1] - b[1], b[0] - a[0], a[0] * b[1] - a[1] * b[0];
        q << a[3] - b[3], b[2] - a[2], a[2] * b[3] - a[3] * b[2];

        if ((p[0]*c[0] + p[1]*c[1] + p[2]) * (q[0]*c[2] + q[1]*c[3] + q[2]) < 0)
            return false;
        if ((p[0]*d[0] + p[1]*d[1] + p[2]) * (q[0]*d[2] + q[1]*d[3] + q[2]) < 0)
            return false;

        p << c[1] - d[1], d[0] - c[0], c[0] * d[1] - c[1] * d[0];
        q << c[3] - d[3], d[2] - c[2], c[2] * d[3] - c[3] * d[2];

        if ((p[0]*a[0] + p[1]*a[1] + p[2]) * (q[0]*a[2] + q[1]*a[3] + q[2]) < 0)
            return false;
        if ((p[0]*b[0] + p[1]*b[1] + p[2]) * (q[0]*b[2] + q[1]*b[3] + q[2]) < 0)
            return false;

        return true;
    }
};

// ---------------------------------------------------------------------------
// MAGSAC estimator extension — adds sigma-consensus constexpr methods
// ---------------------------------------------------------------------------
template <class Solver,
          class _MinimalSolverEngine,
          class _NonMinimalSolverEngine>
class HomLibMAGSACEstimator
    : public HomLibGCEstimator<Solver, _MinimalSolverEngine, _NonMinimalSolverEngine> {
    using Base = HomLibGCEstimator<Solver, _MinimalSolverEngine, _NonMinimalSolverEngine>;

public:
    using Base::residual;

    HomLibMAGSACEstimator(const Solver& solver,
                          const NonMinimalSettings& nm = {},
                          const NonMinimalEstimateFn& nm_solver_fn = {})
        : Base(solver, nm, nm_solver_fn) {}

    inline double residualForScoring(
        const cv::Mat& point_,
        const gcransac::Model& model_) const
    {
        return residual(point_, model_);
    }

    static constexpr double getSigmaQuantile() { return 3.64; }
    static constexpr size_t getDegreesOfFreedom() { return 4; }
    static constexpr double getC() { return 0.25; }
    static constexpr double getGammaFunction() { return 1.0; }
    static constexpr double getUpperIncompleteGammaOfK() { return 0.0036572608340910764; }
    static constexpr double getLowerIncompleteGammaOfK() { return 1.3012265540498875; }
    static constexpr double getChiSquareParamCp() { return 1.0 / (4.0 * getGammaFunction()); }
    static constexpr bool doesNormalizationForNonMinimalFitting() { return true; }
};

} // namespace HomLib

#endif // HOMLIB_HAS_GCRANSAC
