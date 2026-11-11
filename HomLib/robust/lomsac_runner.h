#pragma once

#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <RansacLib/ransac.h>

#include <posedata.hpp>
#include <ransac_estimator.h>
#include <affine_ransac_estimator.h>
#include <orientation_ransac_estimator.h>

namespace HomLib {

template <typename Solver>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac_point(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const NonMinimalSettings& nm,
    const ransac_lib::LORansacOptions& options)
{
    RansacEstimator<Solver> estimator(x, y, solver);
    estimator.nm_settings = nm;

    ransac_lib::RansacStatistics stats;
    PoseData best_model;
    best_model.homography = Eigen::Matrix3d::Identity();
    best_model.focal_length = 0.0;
    best_model.distortion_parameter = 0.0;
    best_model.distortion_parameter2 = 0.0;

    ransac_lib::LocallyOptimizedMSAC<
        PoseData,
        std::vector<PoseData>,
        RansacEstimator<Solver>> lomsac;
    lomsac.EstimateModel(options, estimator, &best_model, &stats);

    return {std::move(best_model), std::move(stats)};
}

template <typename Solver>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac_affine(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A,
    const NonMinimalSettings& nm,
    const ransac_lib::LORansacOptions& options)
{
    AffineRansacEstimator<Solver> estimator(x, y, A, solver, nm);

    ransac_lib::RansacStatistics stats;
    PoseData best_model;
    best_model.homography = Eigen::Matrix3d::Identity();
    best_model.focal_length = 0.0;
    best_model.distortion_parameter = 0.0;
    best_model.distortion_parameter2 = 0.0;

    ransac_lib::LocallyOptimizedMSAC<
        PoseData,
        std::vector<PoseData>,
        AffineRansacEstimator<Solver>> lomsac;
    lomsac.EstimateModel(options, estimator, &best_model, &stats);

    return {std::move(best_model), std::move(stats)};
}

template <typename Solver>
std::pair<PoseData, ransac_lib::RansacStatistics> run_lomsac_orientation(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori,
    const NonMinimalSettings& nm,
    const ransac_lib::LORansacOptions& options)
{
    OrientationRansacEstimator<Solver> estimator(x, y, ori, solver, nm);

    ransac_lib::RansacStatistics stats;
    PoseData best_model;
    best_model.homography = Eigen::Matrix3d::Identity();
    best_model.focal_length = 0.0;
    best_model.distortion_parameter = 0.0;
    best_model.distortion_parameter2 = 0.0;

    ransac_lib::LocallyOptimizedMSAC<
        PoseData,
        std::vector<PoseData>,
        OrientationRansacEstimator<Solver>> lomsac;
    lomsac.EstimateModel(options, estimator, &best_model, &stats);

    return {std::move(best_model), std::move(stats)};
}

} // namespace HomLib
