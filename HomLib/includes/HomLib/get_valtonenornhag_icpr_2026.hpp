// Copyright (c) 2020 Marcus Valtonen Örnhag
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef INCLUDES_HOMLIB_GET_VALTONENORNHAG_ICPR_2026_HPP_
#define INCLUDES_HOMLIB_GET_VALTONENORNHAG_ICPR_2026_HPP_

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "affine_pose_estimator.h"
#include "orientation_pose_estimator.h"
#include "refinement.hpp"

namespace HomLib {
namespace ValtonenOrnhagICPR2026 {
std::vector<HomLib::PoseData> get_affine(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Matrix2d> &A,
        bool extra_check
);
std::vector<HomLib::PoseData> get_ori(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Vector2d> &ori
);

class AffineSolverSingleSided : public AffinePoseEstimator<AffineSolverSingleSided> {
    public:
        static constexpr int minimal_sample_size_value = 2;
        static constexpr bool supports_overdetermined = true;
        static constexpr DistortionCase DISTORTION_CASE = ONE_SIDED_RIGHT;
        bool extra_check = false;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::ValtonenOrnhagICPR2026::get_affine(x, y, A, extra_check);
            HomLib::filter_distortion_range(&output, range);   // a single solution at the minimal sample size
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 2;
        }
        // xu: xd undistorted, yd: yu distorted; false where there is none (see radial.hpp)
        inline bool undistort(const HomLib::PoseData &pose, const Eigen::Vector2d &xd, Eigen::Vector2d &xu) const {
            return HomLib::try_radialundistort(xd, pose.distortion_parameter, xu);
        }
        inline bool distort(const HomLib::PoseData &, const Eigen::Vector2d &yu, Eigen::Vector2d &yd) const {
            yd = yu;
            return true;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided_right(x, y, pose);
        }
    };

class OrientationSolverSingleSided : public OrientationPoseEstimator<OrientationSolverSingleSided> {
    public:
        static constexpr int minimal_sample_size_value = 4;
        static constexpr bool supports_overdetermined = true;
        static constexpr DistortionCase DISTORTION_CASE = ONE_SIDED_RIGHT;
        OrientationSolverSingleSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Vector2d> &ori, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::ValtonenOrnhagICPR2026::get_ori(x, y, ori);
            HomLib::filter_distortion_range(&output, range);   // a single solution at the minimal sample size
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 4;  // Degenerates ?
        }
        // xu: xd undistorted, yd: yu distorted; false where there is none (see radial.hpp)
        inline bool undistort(const HomLib::PoseData &pose, const Eigen::Vector2d &xd, Eigen::Vector2d &xu) const {
            return HomLib::try_radialundistort(xd, pose.distortion_parameter, xu);
        }
        inline bool distort(const HomLib::PoseData &, const Eigen::Vector2d &yu, Eigen::Vector2d &yd) const {
            yd = yu;
            return true;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided_right(x, y, pose);
        }
    };
}
}

#endif  // INCLUDES_HOMLIB_GET_VALTONENORNHAG_ICPR_2026_HPP_
