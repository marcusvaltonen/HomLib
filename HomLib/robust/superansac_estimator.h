#pragma once
#ifdef HOMLIB_HAS_SUPERANSAC

#include <cmath>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include "estimators/abstract_estimator.h"
#include "estimators/solver_homography_four_point.h"
#include "models/model.h"
#include "utils/types.h"

#include "posedata.hpp"
#include "superansac_model_utils.h"
#include "superansac_solver.h"

namespace HomLib {

// SupeRANSAC estimator for point-based HomLib solvers.
// Implements the Estimator virtual interface required by SupeRansac::run().
template <class Solver>
class SuperansacPointEstimator : public superansac::estimator::Estimator {
public:
    explicit SuperansacPointEstimator(const Solver& solver) : solver_(solver) {
        setMinimalSolver(new HomLibMinimalSolver<Solver>(solver));
        setNonMinimalSolver(new HomLibNonMinimalSolver<Solver>(solver));
    }

    bool isWeightingApplicable() const override { return true; }
    double multError() const override { return 1.0; }
    double logAlpha0(size_t w, size_t h, double scalingFactor = 0.5) const override {
        return std::log(1.0 / (w * h * scalingFactor));
    }
    size_t getDegreesOfFreedom() const override { return 2; }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        std::vector<superansac::models::Model>* models) const override
    {
        return minimalSolver->estimateModel(kData, kSample, sampleSize(), *models);
    }

    FORCE_INLINE bool estimateModelNonminimal(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t& kSampleNumber,
        std::vector<superansac::models::Model>* models,
        const double* kWeights = nullptr) const override
    {
        return nonMinimalSolver->estimateModel(kData, kSample, kSampleNumber, *models, kWeights);
    }

    // Distortion-aware residual: undistort source, project through H, distort target.
    FORCE_INLINE double squaredResidual(const double* point,
        const superansac::models::Model& model) const override
    {
        PoseData pose = model_to_pose(model);
        Eigen::Vector2d x(point[0], point[1]);
        Eigen::Vector2d y(point[2], point[3]);

        Eigen::Vector2d z = solver_.undistort(pose, x);
        Eigen::Vector3d Z = pose.homography * z.homogeneous();
        z = Z.hnormalized();
        z = solver_.distort(pose, z);

        return (y - z).squaredNorm();
    }

    FORCE_INLINE double residual(const double* point,
        const superansac::models::Model& model) const override
    {
        return std::sqrt(squaredResidual(point, model));
    }

    void squaredResiduals(
        const DataMatrix& kData,
        const superansac::models::Model& kModel,
        const size_t kStartRow,
        const size_t kCount,
        double* kOut) const override
    {
        PoseData pose = model_to_pose(kModel);
        const size_t cols = kData.cols();
        const double* row = kData.data() + kStartRow * cols;
        for (size_t i = 0; i < kCount; ++i, row += cols) {
            Eigen::Vector2d x(row[0], row[1]);
            Eigen::Vector2d y(row[2], row[3]);

            Eigen::Vector2d z = solver_.undistort(pose, x);
            Eigen::Vector3d Z = pose.homography * z.homogeneous();
            z = Z.hnormalized();
            z = solver_.distort(pose, z);

            kOut[i] = (y - z).squaredNorm();
        }
    }

    FORCE_INLINE bool isValidModel(superansac::models::Model& model,
        const DataMatrix& kData,
        const size_t* kMinimalSample,
        const double kThreshold,
        bool& modelUpdated) const override
    {
        const double det = std::abs(model.getData().block<3, 3>(0, 0).determinant());
        return det > 1e-4 && det < 10000.0;
    }

private:
    Solver solver_;
};

// SupeRANSAC estimator for affine-covariant features.
template <class Solver>
class SuperansacAffineEstimator : public superansac::estimator::Estimator {
public:
    explicit SuperansacAffineEstimator(const Solver& solver) : solver_(solver) {
        setMinimalSolver(new HomLibAffineMinimalSolver<Solver>(solver));
        setNonMinimalSolver(new HomLibNonMinimalSolver<Solver>(solver));
    }

    bool isWeightingApplicable() const override { return true; }
    double multError() const override { return 1.0; }
    double logAlpha0(size_t w, size_t h, double scalingFactor = 0.5) const override {
        return std::log(1.0 / (w * h * scalingFactor));
    }
    size_t getDegreesOfFreedom() const override { return 2; }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        std::vector<superansac::models::Model>* models) const override
    {
        return minimalSolver->estimateModel(kData, kSample, sampleSize(), *models);
    }

    FORCE_INLINE bool estimateModelNonminimal(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t& kSampleNumber,
        std::vector<superansac::models::Model>* models,
        const double* kWeights = nullptr) const override
    {
        return nonMinimalSolver->estimateModel(kData, kSample, kSampleNumber, *models, kWeights);
    }

    FORCE_INLINE double squaredResidual(const double* point,
        const superansac::models::Model& model) const override
    {
        PoseData pose = model_to_pose(model);
        Eigen::Vector2d x(point[0], point[1]);
        Eigen::Vector2d y(point[2], point[3]);

        Eigen::Vector2d z = solver_.undistort(pose, x);
        Eigen::Vector3d Z = pose.homography * z.homogeneous();
        z = Z.hnormalized();
        z = solver_.distort(pose, z);

        return (y - z).squaredNorm();
    }

    FORCE_INLINE double residual(const double* point,
        const superansac::models::Model& model) const override
    {
        return std::sqrt(squaredResidual(point, model));
    }

    void squaredResiduals(
        const DataMatrix& kData,
        const superansac::models::Model& kModel,
        const size_t kStartRow,
        const size_t kCount,
        double* kOut) const override
    {
        PoseData pose = model_to_pose(kModel);
        const size_t cols = kData.cols();
        const double* row = kData.data() + kStartRow * cols;
        for (size_t i = 0; i < kCount; ++i, row += cols) {
            Eigen::Vector2d x(row[0], row[1]);
            Eigen::Vector2d y(row[2], row[3]);

            Eigen::Vector2d z = solver_.undistort(pose, x);
            Eigen::Vector3d Z = pose.homography * z.homogeneous();
            z = Z.hnormalized();
            z = solver_.distort(pose, z);

            kOut[i] = (y - z).squaredNorm();
        }
    }

    FORCE_INLINE bool isValidModel(superansac::models::Model& model,
        const DataMatrix& kData,
        const size_t* kMinimalSample,
        const double kThreshold,
        bool& modelUpdated) const override
    {
        const double det = std::abs(model.getData().block<3, 3>(0, 0).determinant());
        return det > 1e-4 && det < 10000.0;
    }

private:
    Solver solver_;
};

// SupeRANSAC estimator for orientation-covariant features.
template <class Solver>
class SuperansacOrientationEstimator : public superansac::estimator::Estimator {
public:
    explicit SuperansacOrientationEstimator(const Solver& solver) : solver_(solver) {
        setMinimalSolver(new HomLibOrientationMinimalSolver<Solver>(solver));
        setNonMinimalSolver(new HomLibNonMinimalSolver<Solver>(solver));
    }

    bool isWeightingApplicable() const override { return true; }
    double multError() const override { return 1.0; }
    double logAlpha0(size_t w, size_t h, double scalingFactor = 0.5) const override {
        return std::log(1.0 / (w * h * scalingFactor));
    }
    size_t getDegreesOfFreedom() const override { return 2; }

    FORCE_INLINE bool estimateModel(
        const DataMatrix& kData,
        const size_t* kSample,
        std::vector<superansac::models::Model>* models) const override
    {
        return minimalSolver->estimateModel(kData, kSample, sampleSize(), *models);
    }

    FORCE_INLINE bool estimateModelNonminimal(
        const DataMatrix& kData,
        const size_t* kSample,
        const size_t& kSampleNumber,
        std::vector<superansac::models::Model>* models,
        const double* kWeights = nullptr) const override
    {
        return nonMinimalSolver->estimateModel(kData, kSample, kSampleNumber, *models, kWeights);
    }

    FORCE_INLINE double squaredResidual(const double* point,
        const superansac::models::Model& model) const override
    {
        PoseData pose = model_to_pose(model);
        Eigen::Vector2d x(point[0], point[1]);
        Eigen::Vector2d y(point[2], point[3]);

        Eigen::Vector2d z = solver_.undistort(pose, x);
        Eigen::Vector3d Z = pose.homography * z.homogeneous();
        z = Z.hnormalized();
        z = solver_.distort(pose, z);

        return (y - z).squaredNorm();
    }

    FORCE_INLINE double residual(const double* point,
        const superansac::models::Model& model) const override
    {
        return std::sqrt(squaredResidual(point, model));
    }

    void squaredResiduals(
        const DataMatrix& kData,
        const superansac::models::Model& kModel,
        const size_t kStartRow,
        const size_t kCount,
        double* kOut) const override
    {
        PoseData pose = model_to_pose(kModel);
        const size_t cols = kData.cols();
        const double* row = kData.data() + kStartRow * cols;
        for (size_t i = 0; i < kCount; ++i, row += cols) {
            Eigen::Vector2d x(row[0], row[1]);
            Eigen::Vector2d y(row[2], row[3]);

            Eigen::Vector2d z = solver_.undistort(pose, x);
            Eigen::Vector3d Z = pose.homography * z.homogeneous();
            z = Z.hnormalized();
            z = solver_.distort(pose, z);

            kOut[i] = (y - z).squaredNorm();
        }
    }

    FORCE_INLINE bool isValidModel(superansac::models::Model& model,
        const DataMatrix& kData,
        const size_t* kMinimalSample,
        const double kThreshold,
        bool& modelUpdated) const override
    {
        const double det = std::abs(model.getData().block<3, 3>(0, 0).determinant());
        return det > 1e-4 && det < 10000.0;
    }

private:
    Solver solver_;
};

} // namespace HomLib

#endif // HOMLIB_HAS_SUPERANSAC
