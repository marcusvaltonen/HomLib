#pragma once
#ifdef HOMLIB_HAS_SUPERANSAC

#include <vector>
#include <Eigen/Dense>

#include "estimators/abstract_solver.h"
#include "models/model.h"
#include "utils/types.h"
#include "posedata.hpp"
#include "superansac_model_utils.h"

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

// Non-minimal solver: uses SupeRANSAC's DLT HomographyFourPointSolver for the
// initial H estimate, then refines with HomLib's distortion-aware LM.
template <class Solver>
class HomLibNonMinimalSolver : public superansac::estimator::solver::AbstractSolver {
public:
    explicit HomLibNonMinimalSolver(const Solver& solver) : solver_(solver) {}

    bool returnMultipleModels() const override { return false; }
    size_t maximumSolutions() const override { return 1; }
    size_t sampleSize() const override { return 4; }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t kSampleNumber,
        std::vector<superansac::models::Model>& models,
        const double* kWeights = nullptr) const override
    {
        // Use the DLT solver for a good initial H
        if (!dlt_.estimateModel(kData, kSample, kSampleNumber, models, kWeights))
            return false;
        if (models.empty())
            return false;

        // Extract points for HomLib refinement
        std::vector<Eigen::Vector2d> x(kSampleNumber), y(kSampleNumber);
        for (size_t i = 0; i < kSampleNumber; ++i) {
            const size_t idx = kSample ? kSample[i] : i;
            x[i] = Eigen::Vector2d(kData(idx, 0), kData(idx, 1));
            y[i] = Eigen::Vector2d(kData(idx, 2), kData(idx, 3));
        }

        // Convert the DLT result to PoseData, refine, and pack back
        auto& model = models[0];
        PoseData pose;
        pose.homography = model.getData().block<3, 3>(0, 0);
        pose.focal_length = 0.0;
        pose.distortion_parameter = 0.0;
        pose.distortion_parameter2 = 0.0;

        solver_.refine(pose, x, y);

        model = pose_to_model(pose);
        return true;
    }

private:
    Solver solver_;
    superansac::estimator::solver::HomographyFourPointSolver dlt_;
};

} // namespace HomLib

#endif // HOMLIB_HAS_SUPERANSAC
