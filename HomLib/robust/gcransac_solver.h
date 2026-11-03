#pragma once
#ifdef HOMLIB_HAS_GCRANSAC

#include <vector>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "estimators/solver_engine.h"
#include "model.h"
#include "posedata.hpp"
#include "gcransac_model_utils.h"
#include "solver_utils.h"

namespace HomLib {

// ---------------------------------------------------------------------------
// Minimal solver engines — wrap HomLib solvers for GCRANSAC/MAGSAC
// ---------------------------------------------------------------------------

// Point-based minimal solver: reads (x1,y1,x2,y2) from cv::Mat cols 0-3.
template <class Solver, size_t SampleSize>
class GCPointMinimalSolver : public gcransac::estimator::solver::SolverEngine {
public:
    explicit GCPointMinimalSolver(const Solver& solver) : solver_(solver) {}

    static constexpr size_t sampleSize() { return SampleSize; }
    static constexpr size_t maximumSolutions() { return 16; }
    static constexpr bool returnMultipleModels() { return true; }

    OLGA_INLINE bool estimateModel(
        const cv::Mat& data_, const size_t* sample_, size_t sample_number_,
        std::vector<gcransac::Model>& models_,
        const double* weights_ = nullptr) const override
    {
        std::vector<Eigen::Vector2d> x(sample_number_), y(sample_number_);
        const double* ptr = reinterpret_cast<const double*>(data_.data);
        const size_t cols = data_.cols;

        for (size_t i = 0; i < sample_number_; ++i) {
            const double* row = ptr + cols * sample_[i];
            x[i] = Eigen::Vector2d(row[0], row[1]);
            y[i] = Eigen::Vector2d(row[2], row[3]);
        }

        std::vector<PoseData> poses;
        solver_.estimate(x, y, &poses);
        for (const auto& p : poses)
            models_.push_back(pose_to_gc_model(p));
        return !models_.empty();
    }

private:
    Solver solver_;
};

// Affine minimal solver: reads affine params from cols 4-7.
template <class Solver, size_t SampleSize>
class GCAffineMinimalSolver : public gcransac::estimator::solver::SolverEngine {
public:
    explicit GCAffineMinimalSolver(const Solver& solver) : solver_(solver) {}

    static constexpr size_t sampleSize() { return SampleSize; }
    static constexpr size_t maximumSolutions() { return 16; }
    static constexpr bool returnMultipleModels() { return true; }

    OLGA_INLINE bool estimateModel(
        const cv::Mat& data_, const size_t* sample_, size_t sample_number_,
        std::vector<gcransac::Model>& models_,
        const double* weights_ = nullptr) const override
    {
        std::vector<Eigen::Vector2d> x(sample_number_), y(sample_number_);
        std::vector<Eigen::Matrix2d> A(sample_number_);
        const double* ptr = reinterpret_cast<const double*>(data_.data);
        const size_t cols = data_.cols;

        for (size_t i = 0; i < sample_number_; ++i) {
            const double* row = ptr + cols * sample_[i];
            x[i] = Eigen::Vector2d(row[0], row[1]);
            y[i] = Eigen::Vector2d(row[2], row[3]);
            A[i] << row[4], row[5], row[6], row[7];
        }

        std::vector<PoseData> poses;
        solver_.estimate(x, y, A, &poses);
        for (const auto& p : poses)
            models_.push_back(pose_to_gc_model(p));
        return !models_.empty();
    }

private:
    Solver solver_;
};

// Orientation minimal solver: reads orientations from cols 4-5.
template <class Solver, size_t SampleSize>
class GCOrientationMinimalSolver : public gcransac::estimator::solver::SolverEngine {
public:
    explicit GCOrientationMinimalSolver(const Solver& solver) : solver_(solver) {}

    static constexpr size_t sampleSize() { return SampleSize; }
    static constexpr size_t maximumSolutions() { return 16; }
    static constexpr bool returnMultipleModels() { return true; }

    OLGA_INLINE bool estimateModel(
        const cv::Mat& data_, const size_t* sample_, size_t sample_number_,
        std::vector<gcransac::Model>& models_,
        const double* weights_ = nullptr) const override
    {
        std::vector<Eigen::Vector2d> x(sample_number_), y(sample_number_);
        std::vector<Eigen::Vector2d> ori(sample_number_);
        const double* ptr = reinterpret_cast<const double*>(data_.data);
        const size_t cols = data_.cols;

        for (size_t i = 0; i < sample_number_; ++i) {
            const double* row = ptr + cols * sample_[i];
            x[i] = Eigen::Vector2d(row[0], row[1]);
            y[i] = Eigen::Vector2d(row[2], row[3]);
            ori[i] = Eigen::Vector2d(row[4], row[5]);
        }

        std::vector<PoseData> poses;
        solver_.estimate(x, y, ori, &poses);
        for (const auto& p : poses)
            models_.push_back(pose_to_gc_model(p));
        return !models_.empty();
    }

private:
    Solver solver_;
};

// ---------------------------------------------------------------------------
// Non-minimal solver engines — use shared select_and_refine
// ---------------------------------------------------------------------------

// Point-based non-minimal solver.
template <class Solver, size_t SampleSize>
class GCPointNonMinimalSolver : public gcransac::estimator::solver::SolverEngine {
public:
    explicit GCPointNonMinimalSolver(const Solver& solver,
        NonMinimalStrategy strategy = NonMinimalStrategy::Lomsac)
        : solver_(solver), strategy_(strategy) {}

    static constexpr size_t sampleSize() { return SampleSize; }
    static constexpr size_t maximumSolutions() { return 1; }
    static constexpr bool returnMultipleModels() { return false; }

    OLGA_INLINE bool estimateModel(
        const cv::Mat& data_, const size_t* sample_, size_t sample_number_,
        std::vector<gcransac::Model>& models_,
        const double* weights_ = nullptr) const override
    {
        if (sample_number_ < SampleSize)
            return false;

        std::vector<Eigen::Vector2d> x(sample_number_), y(sample_number_);
        const double* ptr = reinterpret_cast<const double*>(data_.data);
        const size_t cols = data_.cols;

        for (size_t i = 0; i < sample_number_; ++i) {
            const size_t idx = sample_ ? sample_[i] : i;
            const double* row = ptr + cols * idx;
            x[i] = Eigen::Vector2d(row[0], row[1]);
            y[i] = Eigen::Vector2d(row[2], row[3]);
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
        if (!detail::non_minimal_dispatch_point(strategy_, solver_, x, y, est_fn, best))
            return false;

        models_.push_back(pose_to_gc_model(best));
        return true;
    }

private:
    Solver solver_;
    NonMinimalStrategy strategy_;
};

// Affine non-minimal solver.
template <class Solver, size_t SampleSize>
class GCAffineNonMinimalSolver : public gcransac::estimator::solver::SolverEngine {
public:
    explicit GCAffineNonMinimalSolver(const Solver& solver,
        NonMinimalStrategy strategy = NonMinimalStrategy::Lomsac)
        : solver_(solver), strategy_(strategy) {}

    static constexpr size_t sampleSize() { return SampleSize; }
    static constexpr size_t maximumSolutions() { return 1; }
    static constexpr bool returnMultipleModels() { return false; }

    OLGA_INLINE bool estimateModel(
        const cv::Mat& data_, const size_t* sample_, size_t sample_number_,
        std::vector<gcransac::Model>& models_,
        const double* weights_ = nullptr) const override
    {
        if (sample_number_ < SampleSize)
            return false;

        std::vector<Eigen::Vector2d> x(sample_number_), y(sample_number_);
        std::vector<Eigen::Matrix2d> A(sample_number_);
        const double* ptr = reinterpret_cast<const double*>(data_.data);
        const size_t cols = data_.cols;

        for (size_t i = 0; i < sample_number_; ++i) {
            const size_t idx = sample_ ? sample_[i] : i;
            const double* row = ptr + cols * idx;
            x[i] = Eigen::Vector2d(row[0], row[1]);
            y[i] = Eigen::Vector2d(row[2], row[3]);
            A[i] << row[4], row[5], row[6], row[7];
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

        models_.push_back(pose_to_gc_model(best));
        return true;
    }

private:
    Solver solver_;
    NonMinimalStrategy strategy_;
};

// Orientation non-minimal solver.
template <class Solver, size_t SampleSize>
class GCOrientationNonMinimalSolver : public gcransac::estimator::solver::SolverEngine {
public:
    explicit GCOrientationNonMinimalSolver(const Solver& solver,
        NonMinimalStrategy strategy = NonMinimalStrategy::Lomsac)
        : solver_(solver), strategy_(strategy) {}

    static constexpr size_t sampleSize() { return SampleSize; }
    static constexpr size_t maximumSolutions() { return 1; }
    static constexpr bool returnMultipleModels() { return false; }

    OLGA_INLINE bool estimateModel(
        const cv::Mat& data_, const size_t* sample_, size_t sample_number_,
        std::vector<gcransac::Model>& models_,
        const double* weights_ = nullptr) const override
    {
        if (sample_number_ < SampleSize)
            return false;

        std::vector<Eigen::Vector2d> x(sample_number_), y(sample_number_);
        std::vector<Eigen::Vector2d> ori(sample_number_);
        const double* ptr = reinterpret_cast<const double*>(data_.data);
        const size_t cols = data_.cols;

        for (size_t i = 0; i < sample_number_; ++i) {
            const size_t idx = sample_ ? sample_[i] : i;
            const double* row = ptr + cols * idx;
            x[i] = Eigen::Vector2d(row[0], row[1]);
            y[i] = Eigen::Vector2d(row[2], row[3]);
            ori[i] = Eigen::Vector2d(row[4], row[5]);
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

        models_.push_back(pose_to_gc_model(best));
        return true;
    }

private:
    Solver solver_;
    NonMinimalStrategy strategy_;
};

} // namespace HomLib

#endif // HOMLIB_HAS_GCRANSAC
