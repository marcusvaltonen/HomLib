#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "distortion_range.h"
#include "solver_id.h"
#include "pose_estimator.h"
#include "affine_pose_estimator.h"
#include "orientation_pose_estimator.h"
#include "refinement.hpp"

namespace HomLib {

namespace detail {

constexpr DistortionCase flipped_case(DistortionCase c) {
    return c == ONE_SIDED_LEFT ? ONE_SIDED_RIGHT : ONE_SIDED_LEFT;
}

inline void flip_poses(std::vector<PoseData>& poses) {
    auto it = std::remove_if(poses.begin(), poses.end(),
        [](const PoseData& p) { return std::abs(p.homography.determinant()) < 1e-15; });
    poses.erase(it, poses.end());
    for (auto& p : poses) {
        p.homography = p.homography.inverse();
        std::swap(p.distortion_parameter, p.distortion_parameter2);
    }
}

template <DistortionCase DC>
inline Eigen::Vector2d flipped_undistort(const PoseData& pose, const Eigen::Vector2d& xd) {
    if constexpr (DC == ONE_SIDED_RIGHT)
        return radialundistort(xd, pose.distortion_parameter);
    else
        return xd;
}

template <DistortionCase DC>
inline Eigen::Vector2d flipped_distort(const PoseData& pose, const Eigen::Vector2d& yu) {
    if constexpr (DC == ONE_SIDED_LEFT)
        return radialdistort(yu, pose.distortion_parameter2);
    else
        return yu;
}

template <DistortionCase DC>
inline void flipped_refine(PoseData& pose, const std::vector<Eigen::Vector2d>& x,
                           const std::vector<Eigen::Vector2d>& y) {
    if constexpr (DC == ONE_SIDED_RIGHT)
        refinement_onesided_right(x, y, pose);
    else
        refinement_onesided(x, y, pose);
}

} // namespace detail

// ---------------------------------------------------------------------------
// Point-based flipped solver
// ---------------------------------------------------------------------------

template <typename InnerSolver>
class FlippedOneSidedSolver : public PoseEstimator<FlippedOneSidedSolver<InnerSolver>> {
    static_assert(InnerSolver::DISTORTION_CASE == ONE_SIDED_LEFT ||
                  InnerSolver::DISTORTION_CASE == ONE_SIDED_RIGHT,
                  "FlippedOneSidedSolver requires a one-sided solver");
public:
    static constexpr int minimal_sample_size_value = InnerSolver::minimal_sample_size_value;
    static constexpr DistortionCase DISTORTION_CASE = detail::flipped_case(InnerSolver::DISTORTION_CASE);

    int solve(const std::vector<Eigen::Vector2d>& x, const std::vector<Eigen::Vector2d>& y,
              std::vector<PoseData>* poses,
              const DistortionRange& range = DistortionRange()) const {
        std::vector<PoseData> output;
        inner_.solve(y, x, &output, range);
        detail::flip_poses(output);
        for (auto& p : output)
            poses->push_back(std::move(p));
        return static_cast<int>(output.size());
    }

    int minimal_sample_size() const { return minimal_sample_size_value; }

    inline Eigen::Vector2d undistort(const PoseData pose, const Eigen::Vector2d& xd) const {
        return detail::flipped_undistort<DISTORTION_CASE>(pose, xd);
    }
    inline Eigen::Vector2d distort(const PoseData pose, const Eigen::Vector2d& yu) const {
        return detail::flipped_distort<DISTORTION_CASE>(pose, yu);
    }
    inline void refine(PoseData& pose, const std::vector<Eigen::Vector2d>& x,
                       const std::vector<Eigen::Vector2d>& y) const {
        detail::flipped_refine<DISTORTION_CASE>(pose, x, y);
    }

    InnerSolver inner_;
};

// ---------------------------------------------------------------------------
// Affine-based flipped solver
// ---------------------------------------------------------------------------

template <typename InnerSolver>
class FlippedAffineOneSidedSolver : public AffinePoseEstimator<FlippedAffineOneSidedSolver<InnerSolver>> {
    static_assert(InnerSolver::DISTORTION_CASE == ONE_SIDED_LEFT ||
                  InnerSolver::DISTORTION_CASE == ONE_SIDED_RIGHT,
                  "FlippedAffineOneSidedSolver requires a one-sided solver");
public:
    static constexpr int minimal_sample_size_value = InnerSolver::minimal_sample_size_value;
    static constexpr DistortionCase DISTORTION_CASE = detail::flipped_case(InnerSolver::DISTORTION_CASE);

    int solve(const std::vector<Eigen::Vector2d>& x, const std::vector<Eigen::Vector2d>& y,
              const std::vector<Eigen::Matrix2d>& A, std::vector<PoseData>* poses,
              const DistortionRange& range = DistortionRange()) const {
        std::vector<Eigen::Matrix2d> A_inv(A.size());
        for (size_t i = 0; i < A.size(); ++i)
            A_inv[i] = A[i].inverse();
        std::vector<PoseData> output;
        inner_.solve(y, x, A_inv, &output, range);
        detail::flip_poses(output);
        for (auto& p : output)
            poses->push_back(std::move(p));
        return static_cast<int>(output.size());
    }

    int minimal_sample_size() const { return minimal_sample_size_value; }

    inline Eigen::Vector2d undistort(const PoseData pose, const Eigen::Vector2d& xd) const {
        return detail::flipped_undistort<DISTORTION_CASE>(pose, xd);
    }
    inline Eigen::Vector2d distort(const PoseData pose, const Eigen::Vector2d& yu) const {
        return detail::flipped_distort<DISTORTION_CASE>(pose, yu);
    }
    inline void refine(PoseData& pose, const std::vector<Eigen::Vector2d>& x,
                       const std::vector<Eigen::Vector2d>& y) const {
        detail::flipped_refine<DISTORTION_CASE>(pose, x, y);
    }

    InnerSolver inner_;
};

// ---------------------------------------------------------------------------
// Orientation-based flipped solver
// ---------------------------------------------------------------------------

template <typename InnerSolver>
class FlippedOrientationOneSidedSolver : public OrientationPoseEstimator<FlippedOrientationOneSidedSolver<InnerSolver>> {
    static_assert(InnerSolver::DISTORTION_CASE == ONE_SIDED_LEFT ||
                  InnerSolver::DISTORTION_CASE == ONE_SIDED_RIGHT,
                  "FlippedOrientationOneSidedSolver requires a one-sided solver");
public:
    static constexpr int minimal_sample_size_value = InnerSolver::minimal_sample_size_value;
    static constexpr DistortionCase DISTORTION_CASE = detail::flipped_case(InnerSolver::DISTORTION_CASE);

    int solve(const std::vector<Eigen::Vector2d>& x, const std::vector<Eigen::Vector2d>& y,
              const std::vector<Eigen::Vector2d>& ori, std::vector<PoseData>* poses,
              const DistortionRange& range = DistortionRange()) const {
        std::vector<Eigen::Vector2d> ori_swapped(ori.size());
        for (size_t i = 0; i < ori.size(); ++i)
            ori_swapped[i] = Eigen::Vector2d(ori[i](1), ori[i](0));
        std::vector<PoseData> output;
        inner_.solve(y, x, ori_swapped, &output, range);
        detail::flip_poses(output);
        for (auto& p : output)
            poses->push_back(std::move(p));
        return static_cast<int>(output.size());
    }

    int minimal_sample_size() const { return minimal_sample_size_value; }

    inline Eigen::Vector2d undistort(const PoseData pose, const Eigen::Vector2d& xd) const {
        return detail::flipped_undistort<DISTORTION_CASE>(pose, xd);
    }
    inline Eigen::Vector2d distort(const PoseData pose, const Eigen::Vector2d& yu) const {
        return detail::flipped_distort<DISTORTION_CASE>(pose, yu);
    }
    inline void refine(PoseData& pose, const std::vector<Eigen::Vector2d>& x,
                       const std::vector<Eigen::Vector2d>& y) const {
        detail::flipped_refine<DISTORTION_CASE>(pose, x, y);
    }

    InnerSolver inner_;
};

} // namespace HomLib
