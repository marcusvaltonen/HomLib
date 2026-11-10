#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <numeric>
#include <Eigen/Dense>
#include "posedata.hpp"

namespace HomLib {

enum class NonMinimalStrategy { Sample, Overdetermined };

namespace detail {

template <class Solver>
inline double compute_squared_residual(
    const Solver& solver, const PoseData& pose,
    const Eigen::Vector2d& x, const Eigen::Vector2d& y)
{
    // max() (an outlier for every backend) where x can't be undistorted, H z is at infinity, H z
    // can't be distorted, or r2 overflows; tested before the operations, so that nothing traps with
    // floating-point exceptions enabled (for a finite H). A NaN would fail every comparison. Not
    // inf: sqrt(r2), 0 * r2 and r2 - r2 stay finite.
    constexpr double kOutlier = std::numeric_limits<double>::max();
    Eigen::Vector2d z, w;
    if (!solver.undistort(pose, x, z))
        return kOutlier;
    const Eigen::Vector3d Z = pose.homography * z.homogeneous();
    if (Z(2) == 0.0 || !solver.distort(pose, Z.hnormalized(), w))
        return kOutlier;
    const double r2 = (y - w).squaredNorm();
    // r2 >= 0: false only for an overflow (inf) or a NaN (a NaN in H); a quiet comparison
    return std::islessequal(r2, kOutlier) ? r2 : kOutlier;
}

template <class Solver>
inline int find_best_pose(
    const Solver& solver,
    const std::vector<PoseData>& poses,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y)
{
    double best_score = std::numeric_limits<double>::max();
    int best_idx = -1;
    for (size_t i = 0; i < poses.size(); ++i) {
        double score = 0.0;
        for (size_t j = 0; j < x.size(); ++j)
            score += compute_squared_residual(solver, poses[i], x[j], y[j]);
        if (score < best_score) {
            best_score = score;
            best_idx = static_cast<int>(i);
        }
    }
    return best_idx;
}

template <class Solver>
inline bool select_best(
    const Solver& solver,
    std::vector<PoseData>& poses,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    PoseData& best)
{
    if (poses.empty())
        return false;

    int best_idx = find_best_pose(solver, poses, x, y);
    if (best_idx < 0)
        return false;

    best = poses[best_idx];
    return true;
}

template <class Solver>
inline void refine_on_inliers(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<size_t>& inliers,
    PoseData& pose)
{
    if (inliers.empty())
        return;
    std::vector<Eigen::Vector2d> xi(inliers.size()), yi(inliers.size());
    for (size_t i = 0; i < inliers.size(); ++i) {
        xi[i] = x[inliers[i]];
        yi[i] = y[inliers[i]];
    }
    solver.refine(pose, xi, yi);
}

inline bool is_valid_model(const Eigen::Matrix3d& H)
{
    return H.norm() > 1e-15;
}

// LOMSAC full non-minimal: all points → estimate_fn → score → no refine.
// estimate_fn signature: (std::vector<PoseData>*) → void
template <class Solver, class EstimateFunc>
inline bool lomsac_non_minimal_full(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    EstimateFunc estimate_fn,
    PoseData& best)
{
    std::vector<PoseData> poses;
    estimate_fn(&poses);

    int best_idx = find_best_pose(solver, poses, x, y);
    if (best_idx < 0)
        return false;

    best = poses[best_idx];
    return true;
}

// Dispatch for point-based non-minimal solving.
// estimate_fn signature: (const Solver&, const vector<Vector2d>&, const vector<Vector2d>&,
//                         const vector<size_t>&, size_t, vector<PoseData>*)
template <class Solver, class EstimateFunc>
inline bool non_minimal_dispatch_point(
    NonMinimalStrategy strategy,
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    EstimateFunc estimate_fn,
    PoseData& best,
    bool refine_non_minimal = true)
{
    const size_t n = x.size();
    const size_t ss = static_cast<size_t>(solver.minimal_sample_size());
    if (n < ss)
        return false;

    std::vector<size_t> idx(ss);
    for (size_t i = 0; i < ss; ++i)
        idx[i] = i * n / ss;
    std::vector<PoseData> poses;
    estimate_fn(solver, x, y, idx, ss, &poses);
    bool found = select_best(solver, poses, x, y, best);

    if (found && refine_non_minimal)
        solver.refine(best, x, y);
    return found;
}

// Dispatch for full non-minimal solving (affine/orientation — uses all points).
template <class Solver, class EstimateFunc>
inline bool non_minimal_dispatch_full(
    NonMinimalStrategy strategy,
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    EstimateFunc estimate_fn,
    PoseData& best,
    bool refine_non_minimal = true)
{
    const size_t n = x.size();
    if (n < static_cast<size_t>(solver.minimal_sample_size()))
        return false;

    bool found = false;

    switch (strategy) {
        case NonMinimalStrategy::Sample:
        case NonMinimalStrategy::Overdetermined: {
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::vector<PoseData> poses;
            estimate_fn(solver, x, y, idx, n, &poses);
            found = select_best(solver, poses, x, y, best);
            break;
        }
    }

    if (found && refine_non_minimal)
        solver.refine(best, x, y);
    return found;
}

} // namespace detail
} // namespace HomLib
