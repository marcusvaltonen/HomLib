#pragma once

#include <vector>
#include <limits>
#include <random>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include "posedata.hpp"

namespace HomLib {

enum class NonMinimalStrategy { Lomsac, SelectAndRefine, MultiSample };

namespace detail {

template <class Solver>
inline double compute_squared_residual(
    const Solver& solver, const PoseData& pose,
    const Eigen::Vector2d& x, const Eigen::Vector2d& y)
{
    Eigen::Vector2d z = solver.undistort(pose, x);
    Eigen::Vector3d Z = pose.homography * z.homogeneous();
    z = Z.hnormalized();
    z = solver.distort(pose, z);
    return (y - z).squaredNorm();
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
inline bool select_and_refine(
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
    solver.refine(best, x, y);
    return true;
}

// Try multiple random minimal-sample subsets from the inlier set,
// pick the one with the lowest total residual, then refine on all points.
// EstimateFunc: callable(solver, xs, ys, indices, &poses) that calls the
// appropriate solver.estimate() variant with any extra data.
template <class Solver, class EstimateFunc>
inline bool multi_sample_estimate(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    EstimateFunc estimate_fn,
    PoseData& best,
    int n_trials = 10)
{
    const size_t n = x.size();
    const size_t ss = static_cast<size_t>(solver.minimal_sample_size());
    if (n < ss)
        return false;

    thread_local std::mt19937 rng{std::random_device{}()};
    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

    double best_score = std::numeric_limits<double>::max();
    bool found = false;

    for (int t = 0; t < n_trials; ++t) {
        std::shuffle(indices.begin(), indices.end(), rng);

        std::vector<PoseData> poses;
        estimate_fn(solver, x, y, indices, ss, &poses);
        if (poses.empty())
            continue;

        int idx = find_best_pose(solver, poses, x, y);
        if (idx < 0)
            continue;

        double score = 0.0;
        for (size_t j = 0; j < n; ++j)
            score += compute_squared_residual(solver, poses[idx], x[j], y[j]);

        if (score < best_score) {
            best_score = score;
            best = poses[idx];
            found = true;
        }
    }

    if (!found)
        return false;

    solver.refine(best, x, y);
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

// LOMSAC point non-minimal: first min_sample_size points → estimate → score → no refine.
template <class Solver>
inline bool lomsac_non_minimal_point(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    PoseData& best)
{
    const size_t ss = static_cast<size_t>(solver.minimal_sample_size());
    if (x.size() < ss)
        return false;

    std::vector<Eigen::Vector2d> xs(x.begin(), x.begin() + ss);
    std::vector<Eigen::Vector2d> ys(y.begin(), y.begin() + ss);

    std::vector<PoseData> poses;
    solver.estimate(xs, ys, &poses);

    int best_idx = find_best_pose(solver, poses, x, y);
    if (best_idx < 0)
        return false;

    best = poses[best_idx];
    return true;
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
    PoseData& best)
{
    const size_t n = x.size();
    const size_t ss = static_cast<size_t>(solver.minimal_sample_size());
    if (n < ss)
        return false;

    switch (strategy) {
        case NonMinimalStrategy::Lomsac: {
            // Shuffle before taking first-N: backends (GC-RANSAC, SupeRANSAC)
            // may provide inlier indices in fixed order, unlike LOMSAC which
            // provides random samples from RansacLib.
            thread_local std::mt19937 rng{std::random_device{}()};
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng);
            idx.resize(ss);
            std::vector<PoseData> poses;
            estimate_fn(solver, x, y, idx, ss, &poses);
            int best_idx = find_best_pose(solver, poses, x, y);
            if (best_idx < 0)
                return false;
            best = poses[best_idx];
            return true;
        }
        case NonMinimalStrategy::SelectAndRefine: {
            // Estimate on first-N (shuffled) like Lomsac, then refine on all.
            // Point solvers handle small samples better than large overdetermined
            // systems; the Gauss-Newton refine step recovers accuracy.
            thread_local std::mt19937 rng2{std::random_device{}()};
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng2);
            idx.resize(ss);
            std::vector<PoseData> poses;
            estimate_fn(solver, x, y, idx, ss, &poses);
            return select_and_refine(solver, poses, x, y, best);
        }
        case NonMinimalStrategy::MultiSample:
            return multi_sample_estimate(solver, x, y, estimate_fn, best);
    }
    return false;
}

// Dispatch for full non-minimal solving (affine/orientation — uses all points).
template <class Solver, class EstimateFunc>
inline bool non_minimal_dispatch_full(
    NonMinimalStrategy strategy,
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    EstimateFunc estimate_fn,
    PoseData& best)
{
    const size_t n = x.size();
    if (n < static_cast<size_t>(solver.minimal_sample_size()))
        return false;

    switch (strategy) {
        case NonMinimalStrategy::Lomsac: {
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::vector<PoseData> poses;
            estimate_fn(solver, x, y, idx, n, &poses);
            int best_idx = find_best_pose(solver, poses, x, y);
            if (best_idx < 0)
                return false;
            best = poses[best_idx];
            return true;
        }
        case NonMinimalStrategy::SelectAndRefine: {
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::vector<PoseData> poses;
            estimate_fn(solver, x, y, idx, n, &poses);
            return select_and_refine(solver, poses, x, y, best);
        }
        case NonMinimalStrategy::MultiSample:
            return multi_sample_estimate(solver, x, y, estimate_fn, best);
    }
    return false;
}

} // namespace detail
} // namespace HomLib
