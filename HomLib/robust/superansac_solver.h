#pragma once
#ifdef HOMLIB_HAS_SUPERANSAC

#include <limits>
#include <vector>
#include <Eigen/Dense>

#include "estimators/abstract_solver.h"
#include "models/model.h"
#include "utils/types.h"
#include "posedata.hpp"
#include "superansac_model_utils.h"
#include "solver_utils.h"

namespace HomLib {

// Minimal solver: extracts point correspondences from SupeRANSAC's DataMatrix,
// calls HomLib's solver, packs results into 4x4 models.
template <class Solver>
class HomLibMinimalSolver : public superansac::estimator::solver::AbstractSolver {
public:
    explicit HomLibMinimalSolver(const Solver& solver) : solver_(solver) {}

    bool returnMultipleModels() const override { return true; }
    size_t maximumSolutions() const override { return 16; }
    size_t sampleSize() const override { return solver_.minimal_sample_size(); }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t kSampleNumber,
        std::vector<superansac::models::Model>& models,
        const double* kWeights = nullptr) const override
    {
        std::vector<Eigen::Vector2d> x(kSampleNumber), y(kSampleNumber);
        for (size_t i = 0; i < kSampleNumber; ++i) {
            const size_t idx = kSample ? kSample[i] : i;
            x[i] = Eigen::Vector2d(kData(idx, 0), kData(idx, 1));
            y[i] = Eigen::Vector2d(kData(idx, 2), kData(idx, 3));
        }

        std::vector<PoseData> poses;
        solver_.estimate(x, y, &poses);

        for (const auto& pose : poses)
            models.push_back(pose_to_model(pose));
        return !models.empty();
    }

private:
    Solver solver_;
};

// Minimal solver for affine features: also extracts 2x2 affine matrices from cols 4-7.
template <class Solver>
class HomLibAffineMinimalSolver : public superansac::estimator::solver::AbstractSolver {
public:
    explicit HomLibAffineMinimalSolver(const Solver& solver) : solver_(solver) {}

    bool returnMultipleModels() const override { return true; }
    size_t maximumSolutions() const override { return 16; }
    size_t sampleSize() const override { return solver_.minimal_sample_size(); }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t kSampleNumber,
        std::vector<superansac::models::Model>& models,
        const double* kWeights = nullptr) const override
    {
        std::vector<Eigen::Vector2d> x(kSampleNumber), y(kSampleNumber);
        std::vector<Eigen::Matrix2d> A(kSampleNumber);
        for (size_t i = 0; i < kSampleNumber; ++i) {
            const size_t idx = kSample ? kSample[i] : i;
            x[i] = Eigen::Vector2d(kData(idx, 0), kData(idx, 1));
            y[i] = Eigen::Vector2d(kData(idx, 2), kData(idx, 3));
            A[i] << kData(idx, 4), kData(idx, 5),
                     kData(idx, 6), kData(idx, 7);
        }

        std::vector<PoseData> poses;
        solver_.estimate(x, y, A, &poses);

        for (const auto& pose : poses)
            models.push_back(pose_to_model(pose));
        return !models.empty();
    }

private:
    Solver solver_;
};

// Minimal solver for orientation features: extracts orientations from cols 4-5.
template <class Solver>
class HomLibOrientationMinimalSolver : public superansac::estimator::solver::AbstractSolver {
public:
    explicit HomLibOrientationMinimalSolver(const Solver& solver) : solver_(solver) {}

    bool returnMultipleModels() const override { return true; }
    size_t maximumSolutions() const override { return 16; }
    size_t sampleSize() const override { return solver_.minimal_sample_size(); }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t kSampleNumber,
        std::vector<superansac::models::Model>& models,
        const double* kWeights = nullptr) const override
    {
        std::vector<Eigen::Vector2d> x(kSampleNumber), y(kSampleNumber);
        std::vector<Eigen::Vector2d> ori(kSampleNumber);
        for (size_t i = 0; i < kSampleNumber; ++i) {
            const size_t idx = kSample ? kSample[i] : i;
            x[i] = Eigen::Vector2d(kData(idx, 0), kData(idx, 1));
            y[i] = Eigen::Vector2d(kData(idx, 2), kData(idx, 3));
            ori[i] = Eigen::Vector2d(kData(idx, 4), kData(idx, 5));
        }

        std::vector<PoseData> poses;
        solver_.estimate(x, y, ori, &poses);

        for (const auto& pose : poses)
            models.push_back(pose_to_model(pose));
        return !models.empty();
    }

private:
    Solver solver_;
};

// Point non-minimal solver: dispatches via shared non_minimal_dispatch_point.
template <class Solver>
class HomLibPointNonMinimalSolver : public superansac::estimator::solver::AbstractSolver {
public:
    explicit HomLibPointNonMinimalSolver(const Solver& solver,
        NonMinimalStrategy strategy = NonMinimalStrategy::Sample,
        NonMinimalEstimateFn non_minimal_fn = {})
        : solver_(solver), strategy_(strategy), non_minimal_fn_(std::move(non_minimal_fn)) {}

    bool returnMultipleModels() const override { return false; }
    size_t maximumSolutions() const override { return 1; }
    size_t sampleSize() const override { return solver_.minimal_sample_size(); }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t kSampleNumber,
        std::vector<superansac::models::Model>& models,
        const double* kWeights = nullptr) const override
    {
        if (kSampleNumber < static_cast<size_t>(solver_.minimal_sample_size()))
            return false;

        std::vector<Eigen::Vector2d> x(kSampleNumber), y(kSampleNumber);
        for (size_t i = 0; i < kSampleNumber; ++i) {
            const size_t idx = kSample ? kSample[i] : i;
            x[i] = Eigen::Vector2d(kData(idx, 0), kData(idx, 1));
            y[i] = Eigen::Vector2d(kData(idx, 2), kData(idx, 3));
        }

        auto est_fn = [](const Solver& s,
                         const std::vector<Eigen::Vector2d>& xv,
                         const std::vector<Eigen::Vector2d>& yv,
                         const std::vector<size_t>& idx, size_t ss,
                         std::vector<PoseData>* poses) {
            std::vector<Eigen::Vector2d> xs(ss), ys(ss);
            for (size_t i = 0; i < ss; ++i) { xs[i] = xv[idx[i]]; ys[i] = yv[idx[i]]; }
            s.estimate(xs, ys, poses);
        };

        PoseData best;
        if (!detail::non_minimal_dispatch_point(strategy_, solver_, x, y, est_fn, best, -1, non_minimal_fn_))
            return false;
        models.push_back(pose_to_model(best));
        return true;
    }

private:
    Solver solver_;
    NonMinimalStrategy strategy_;
    NonMinimalEstimateFn non_minimal_fn_;
};

// Affine non-minimal solver: dispatches via shared non_minimal_dispatch_full.
template <class Solver>
class HomLibAffineNonMinimalSolver : public superansac::estimator::solver::AbstractSolver {
public:
    explicit HomLibAffineNonMinimalSolver(const Solver& solver,
        NonMinimalStrategy strategy = NonMinimalStrategy::Sample,
        NonMinimalEstimateFn = {})
        : solver_(solver), strategy_(strategy) {}

    bool returnMultipleModels() const override { return false; }
    size_t maximumSolutions() const override { return 1; }
    size_t sampleSize() const override { return solver_.minimal_sample_size(); }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t kSampleNumber,
        std::vector<superansac::models::Model>& models,
        const double* kWeights = nullptr) const override
    {
        if (kSampleNumber < static_cast<size_t>(solver_.minimal_sample_size()))
            return false;

        std::vector<Eigen::Vector2d> x(kSampleNumber), y(kSampleNumber);
        std::vector<Eigen::Matrix2d> A(kSampleNumber);
        for (size_t i = 0; i < kSampleNumber; ++i) {
            const size_t idx = kSample ? kSample[i] : i;
            x[i] = Eigen::Vector2d(kData(idx, 0), kData(idx, 1));
            y[i] = Eigen::Vector2d(kData(idx, 2), kData(idx, 3));
            A[i] << kData(idx, 4), kData(idx, 5),
                     kData(idx, 6), kData(idx, 7);
        }

        auto est_fn = [&A](const Solver& s,
                           const std::vector<Eigen::Vector2d>& xv,
                           const std::vector<Eigen::Vector2d>& yv,
                           const std::vector<size_t>& idx, size_t ss,
                           std::vector<PoseData>* poses) {
            std::vector<Eigen::Vector2d> xs(ss), ys(ss);
            std::vector<Eigen::Matrix2d> As(ss);
            for (size_t i = 0; i < ss; ++i) {
                xs[i] = xv[idx[i]]; ys[i] = yv[idx[i]]; As[i] = A[idx[i]];
            }
            s.estimate(xs, ys, As, poses);
        };

        PoseData best;
        if (!detail::non_minimal_dispatch_full(strategy_, solver_, x, y, est_fn, best))
            return false;
        models.push_back(pose_to_model(best));
        return true;
    }

private:
    Solver solver_;
    NonMinimalStrategy strategy_;
};

// Orientation non-minimal solver: dispatches via shared non_minimal_dispatch_full.
template <class Solver>
class HomLibOrientationNonMinimalSolver : public superansac::estimator::solver::AbstractSolver {
public:
    explicit HomLibOrientationNonMinimalSolver(const Solver& solver,
        NonMinimalStrategy strategy = NonMinimalStrategy::Sample,
        NonMinimalEstimateFn = {})
        : solver_(solver), strategy_(strategy) {}

    bool returnMultipleModels() const override { return false; }
    size_t maximumSolutions() const override { return 1; }
    size_t sampleSize() const override { return solver_.minimal_sample_size(); }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t kSampleNumber,
        std::vector<superansac::models::Model>& models,
        const double* kWeights = nullptr) const override
    {
        if (kSampleNumber < static_cast<size_t>(solver_.minimal_sample_size()))
            return false;

        std::vector<Eigen::Vector2d> x(kSampleNumber), y(kSampleNumber);
        std::vector<Eigen::Vector2d> ori(kSampleNumber);
        for (size_t i = 0; i < kSampleNumber; ++i) {
            const size_t idx = kSample ? kSample[i] : i;
            x[i] = Eigen::Vector2d(kData(idx, 0), kData(idx, 1));
            y[i] = Eigen::Vector2d(kData(idx, 2), kData(idx, 3));
            ori[i] = Eigen::Vector2d(kData(idx, 4), kData(idx, 5));
        }

        auto est_fn = [&ori](const Solver& s,
                             const std::vector<Eigen::Vector2d>& xv,
                             const std::vector<Eigen::Vector2d>& yv,
                             const std::vector<size_t>& idx, size_t ss,
                             std::vector<PoseData>* poses) {
            std::vector<Eigen::Vector2d> xs(ss), ys(ss), os(ss);
            for (size_t i = 0; i < ss; ++i) {
                xs[i] = xv[idx[i]]; ys[i] = yv[idx[i]]; os[i] = ori[idx[i]];
            }
            s.estimate(xs, ys, os, poses);
        };

        PoseData best;
        if (!detail::non_minimal_dispatch_full(strategy_, solver_, x, y, est_fn, best))
            return false;
        models.push_back(pose_to_model(best));
        return true;
    }

private:
    Solver solver_;
    NonMinimalStrategy strategy_;
};

} // namespace HomLib

#endif // HOMLIB_HAS_SUPERANSAC
