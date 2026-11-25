#pragma once

#include <algorithm>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <RansacLib/ransac.h>

#include <posedata.hpp>
#include <ransac_estimator.h>
#include <affine_ransac_estimator.h>
#include <orientation_ransac_estimator.h>
#include <scale_orientation_ransac_estimator.h>

namespace HomLib {

// LO-MSAC on the data of the estimator. RansacLib's final least-squares fit
// (final_least_squares_) also runs when no hypothesis was found: it then refines
// and scores the starting model (the identity), and returns it with the points
// that it happens to fit as inliers. So the fit runs here instead, as RansacLib
// does it, but only on a model that was found.
template <typename Estimator>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac(
    const Estimator& estimator,
    const ransac_lib::LORansacOptions& options)
{
    ransac_lib::LORansacOptions lo_options = options;
    lo_options.final_least_squares_ = false;

    ransac_lib::RansacStatistics stats;
    PoseData best_model;
    best_model.homography = Eigen::Matrix3d::Identity();
    best_model.focal_length = 0.0;
    best_model.distortion_parameter = 0.0;
    best_model.distortion_parameter2 = 0.0;

    ransac_lib::LocallyOptimizedMSAC<
        PoseData,
        std::vector<PoseData>,
        Estimator> lomsac;
    lomsac.EstimateModel(lo_options, estimator, &best_model, &stats);
    if (!options.final_least_squares_ ||
        stats.best_model_score == std::numeric_limits<double>::max())
        return {std::move(best_model), std::move(stats)};

    // Refine on the inliers, and keep the refined model if its MSAC score is lower
    PoseData refined_model = best_model;
    estimator.LeastSquares(stats.inlier_indices, &refined_model);
    const double threshold = options.squared_inlier_threshold_;
    double score = 0.0;
    std::vector<int> inliers;
    for (int i = 0; i < estimator.num_data(); ++i) {
        const double squared_error = estimator.EvaluateModelOnPoint(refined_model, i);
        score += std::min(squared_error, threshold);
        if (squared_error < threshold)
            inliers.push_back(i);
    }
    if (score < stats.best_model_score) {
        best_model = std::move(refined_model);
        stats.best_model_score = score;
        stats.best_num_inliers = static_cast<int>(inliers.size());
        stats.inlier_ratio = static_cast<double>(inliers.size()) / estimator.num_data();
        stats.inlier_indices = std::move(inliers);
    }
    return {std::move(best_model), std::move(stats)};
}

template <typename Solver>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac_point(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const NonMinimalSettings& nm,
    const ransac_lib::LORansacOptions& options,
    const NonMinimalEstimateFn& nm_solver_fn = {})
{
    RansacEstimator<Solver> estimator(x, y, solver);
    estimator.nm_settings = nm;
    estimator.nm_solver_fn = nm_solver_fn;

    return run_lomsac(estimator, options);
}

template <typename Solver>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac_affine(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A,
    const NonMinimalSettings& nm,
    const ransac_lib::LORansacOptions& options,
    const NonMinimalEstimateFn& nm_solver_fn = {})
{
    AffineRansacEstimator<Solver> estimator(x, y, A, solver, nm);
    estimator.nm_solver_fn = nm_solver_fn;

    return run_lomsac(estimator, options);
}

template <typename Solver>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac_orientation(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori,
    const NonMinimalSettings& nm,
    const ransac_lib::LORansacOptions& options,
    const NonMinimalEstimateFn& nm_solver_fn = {})
{
    OrientationRansacEstimator<Solver> estimator(x, y, ori, solver, nm);
    estimator.nm_solver_fn = nm_solver_fn;

    return run_lomsac(estimator, options);
}

template <typename Solver>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac_scale_orientation(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori,
    const std::vector<Eigen::Vector2d>& scale,
    const NonMinimalSettings& nm,
    const ransac_lib::LORansacOptions& options,
    const NonMinimalEstimateFn& nm_solver_fn = {})
{
    ScaleOrientationRansacEstimator<Solver> estimator(x, y, ori, scale, solver, nm);
    estimator.nm_solver_fn = nm_solver_fn;

    ransac_lib::RansacStatistics stats;
    PoseData best_model;
    best_model.homography = Eigen::Matrix3d::Identity();
    best_model.focal_length = 0.0;
    best_model.distortion_parameter = 0.0;
    best_model.distortion_parameter2 = 0.0;

    ransac_lib::LocallyOptimizedMSAC<
        PoseData,
        std::vector<PoseData>,
        ScaleOrientationRansacEstimator<Solver>> lomsac;
    lomsac.EstimateModel(options, estimator, &best_model, &stats);

    return {std::move(best_model), std::move(stats)};
}

} // namespace HomLib
