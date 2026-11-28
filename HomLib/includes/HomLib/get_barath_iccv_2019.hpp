// Copyright (C) 2019 Czech Technical University.
// All rights reserved.

#ifndef INCLUDES_HOMLIB_GET_BARATH_ICCV_2019_HPP_
#define INCLUDES_HOMLIB_GET_BARATH_ICCV_2019_HPP_

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "scale_orientation_pose_estimator.h"
#include "orientation_pose_estimator.h"
#include "refinement.hpp"

namespace HomLib {
namespace BarathICCV2019 {

std::vector<HomLib::PoseData> get_sift(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Vector2d> &ori,
        const std::vector<Eigen::Vector2d> &scale
);

std::vector<HomLib::PoseData> get_ori(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Vector2d> &ori
);

class SIFTSolver : public ScaleOrientationPoseEstimator<SIFTSolver> {
    public:
        static constexpr int minimal_sample_size_value = 2;
        static constexpr bool supports_overdetermined = false;
        static constexpr DistortionCase DISTORTION_CASE = NO_DISTORTION;
        SIFTSolver() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Vector2d> &ori, const std::vector<Eigen::Vector2d> &scale, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange & = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::BarathICCV2019::get_sift(x, y, ori, scale);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 2;
        }
        inline bool undistort(const HomLib::PoseData &, const Eigen::Vector2d &xd, Eigen::Vector2d &xu) const {
            xu = xd;
            return true;
        }
        inline bool distort(const HomLib::PoseData &, const Eigen::Vector2d &yu, Eigen::Vector2d &yd) const {
            yd = yu;
            return true;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_no_dist(x, y, pose);
        }
    };

class OrientationSolver : public OrientationPoseEstimator<OrientationSolver> {
    public:
        static constexpr int minimal_sample_size_value = 3;
        static constexpr bool supports_overdetermined = true;
        static constexpr DistortionCase DISTORTION_CASE = NO_DISTORTION;
        OrientationSolver() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Vector2d> &ori, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange & = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::BarathICCV2019::get_ori(x, y, ori);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 3;
        }
        inline bool undistort(const HomLib::PoseData &, const Eigen::Vector2d &xd, Eigen::Vector2d &xu) const {
            xu = xd;
            return true;
        }
        inline bool distort(const HomLib::PoseData &, const Eigen::Vector2d &yu, Eigen::Vector2d &yd) const {
            yd = yu;
            return true;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_no_dist(x, y, pose);
        }
    };

}  // namespace BarathICCV2019
}  // namespace HomLib

#endif  // INCLUDES_HOMLIB_GET_BARATH_ICCV_2019_HPP_
