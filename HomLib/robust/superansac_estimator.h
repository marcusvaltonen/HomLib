#pragma once
#ifdef HOMLIB_HAS_SUPERANSAC

#include <cmath>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include "estimators/abstract_estimator.h"
#include "models/model.h"
#include "utils/types.h"

#include "posedata.hpp"
#include "superansac_model_utils.h"
#include "superansac_solver.h"
#include "solver_utils.h"

namespace HomLib {

// SupeRANSAC estimator for point-based HomLib solvers.
// Implements the Estimator virtual interface required by SupeRansac::run().
template <class Solver>
class SuperansacPointEstimator : public superansac::estimator::Estimator {
public:
    explicit SuperansacPointEstimator(const Solver& solver,
        const NonMinimalSettings& nm = {},
        const NonMinimalEstimateFn& nm_solver_fn = {})
        : solver_(solver) {
        setMinimalSolver(new HomLibMinimalSolver<Solver>(solver));
        setNonMinimalSolver(new HomLibPointNonMinimalSolver<Solver>(solver, nm, nm_solver_fn));
    }

    // superansac's AbstractSolver has no virtual destructor, so deleting the solvers through
    // unique_ptr<AbstractSolver> is undefined behaviour (clang compiles it to a trap on Apple platforms): delete them here.
    ~SuperansacPointEstimator() override {
        delete static_cast<HomLibMinimalSolver<Solver>*>(minimalSolver.release());
        delete static_cast<HomLibPointNonMinimalSolver<Solver>*>(nonMinimalSolver.release());
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
        return detail::compute_squared_residual(solver_, pose, x, y);
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
            kOut[i] = detail::compute_squared_residual(solver_, pose, x, y);
        }
    }

    FORCE_INLINE bool isValidModel(superansac::models::Model& model,
        const DataMatrix& kData,
        const size_t* kMinimalSample,
        const double kThreshold,
        bool& modelUpdated) const override
    {
        return detail::is_valid_model(model.getData().block<3, 3>(0, 0));
    }

private:
    Solver solver_;
};

// SupeRANSAC estimator for affine-covariant features.
template <class Solver>
class SuperansacAffineEstimator : public superansac::estimator::Estimator {
public:
    explicit SuperansacAffineEstimator(const Solver& solver,
        const NonMinimalSettings& nm = {},
        const NonMinimalEstimateFn& nm_solver_fn = {}) : solver_(solver) {
        setMinimalSolver(new HomLibAffineMinimalSolver<Solver>(solver));
        setNonMinimalSolver(new HomLibAffineNonMinimalSolver<Solver>(solver, nm, nm_solver_fn));
    }

    // superansac's AbstractSolver has no virtual destructor, so deleting the solvers through
    // unique_ptr<AbstractSolver> is undefined behaviour (clang compiles it to a trap on Apple platforms): delete them here.
    ~SuperansacAffineEstimator() override {
        delete static_cast<HomLibAffineMinimalSolver<Solver>*>(minimalSolver.release());
        delete static_cast<HomLibAffineNonMinimalSolver<Solver>*>(nonMinimalSolver.release());
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
        return detail::compute_squared_residual(solver_, pose, x, y);
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
            kOut[i] = detail::compute_squared_residual(solver_, pose, x, y);
        }
    }

    FORCE_INLINE bool isValidModel(superansac::models::Model& model,
        const DataMatrix& kData,
        const size_t* kMinimalSample,
        const double kThreshold,
        bool& modelUpdated) const override
    {
        return detail::is_valid_model(model.getData().block<3, 3>(0, 0));
    }

private:
    Solver solver_;
};

// SupeRANSAC estimator for orientation-covariant features.
template <class Solver>
class SuperansacOrientationEstimator : public superansac::estimator::Estimator {
public:
    explicit SuperansacOrientationEstimator(const Solver& solver,
        const NonMinimalSettings& nm = {},
        const NonMinimalEstimateFn& nm_solver_fn = {}) : solver_(solver) {
        setMinimalSolver(new HomLibOrientationMinimalSolver<Solver>(solver));
        setNonMinimalSolver(new HomLibOrientationNonMinimalSolver<Solver>(solver, nm, nm_solver_fn));
    }

    // superansac's AbstractSolver has no virtual destructor, so deleting the solvers through
    // unique_ptr<AbstractSolver> is undefined behaviour (clang compiles it to a trap on Apple platforms): delete them here.
    ~SuperansacOrientationEstimator() override {
        delete static_cast<HomLibOrientationMinimalSolver<Solver>*>(minimalSolver.release());
        delete static_cast<HomLibOrientationNonMinimalSolver<Solver>*>(nonMinimalSolver.release());
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
        return detail::compute_squared_residual(solver_, pose, x, y);
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
            kOut[i] = detail::compute_squared_residual(solver_, pose, x, y);
        }
    }

    FORCE_INLINE bool isValidModel(superansac::models::Model& model,
        const DataMatrix& kData,
        const size_t* kMinimalSample,
        const double kThreshold,
        bool& modelUpdated) const override
    {
        return detail::is_valid_model(model.getData().block<3, 3>(0, 0));
    }

private:
    Solver solver_;
};

} // namespace HomLib

#endif // HOMLIB_HAS_SUPERANSAC
