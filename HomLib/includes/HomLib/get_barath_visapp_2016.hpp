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

#ifndef INCLUDES_HOMLIB_GET_BARATH_VISAPP_2016_HPP_
#define INCLUDES_HOMLIB_GET_BARATH_VISAPP_2016_HPP_

#include <Eigen/Dense>
#include <limits>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "affine_pose_estimator.h"
#include "orientation_pose_estimator.h"
#include "pose_estimator.h"
#include "refinement.hpp"

namespace HomLib {
namespace BarathVISAPP2016 {
std::vector<HomLib::PoseData> get_affine(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Matrix2d> &A,
        bool use_qr = false
);

std::vector<HomLib::PoseData> get_haf(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Matrix2d> &A,
        const Eigen::Matrix3d &F
);

std::vector<HomLib::PoseData> get_3pt(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const Eigen::Matrix3d &F
);

class AffineSolver : public AffinePoseEstimator<AffineSolver> {
    public:
        static constexpr int minimal_sample_size_value = 2;
        static constexpr bool supports_overdetermined = true;
        static constexpr DistortionCase DISTORTION_CASE = NO_DISTORTION;
        AffineSolver() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange & = HomLib::DistortionRange()) const {   // no distortion (k = 0 passes)
            std::vector<HomLib::PoseData> output = HomLib::BarathVISAPP2016::get_affine(x, y, A);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 2;
        }
        // xu: xd undistorted, yd: yu distorted; false where there is none (see radial.hpp)
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

class HAFSolver : public AffinePoseEstimator<HAFSolver> {
    public:
        static constexpr int minimal_sample_size_value = 1;
        static constexpr bool supports_overdetermined = true;
        static constexpr DistortionCase DISTORTION_CASE = NO_DISTORTION;
        Eigen::Matrix3d F_;
        explicit HAFSolver(const Eigen::Matrix3d &F) : F_(F) {}
        HAFSolver() : F_(Eigen::Matrix3d::Constant(std::numeric_limits<double>::quiet_NaN())) {}
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange & = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::BarathVISAPP2016::get_haf(x, y, A, F_);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        // Override estimate() to transform F_ alongside the point normalization.
        // estimate_normalized() scales points by S = diag(1/f0,1/f0,1); F must follow
        // as F' = S^{-T} F S^{-1} = diag(f0,f0,1) F diag(f0,f0,1).
        int estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_,
                     const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses) const {
            std::vector<Eigen::Vector2d> x = x_, y = y_;
            double f0;
            poses->clear();
            if (!normalize(x, y, &f0))
                return 0;
            Eigen::Matrix3d Fn = Eigen::Vector3d(f0, f0, 1).asDiagonal() * F_ * Eigen::Vector3d(f0, f0, 1).asDiagonal();
            std::vector<HomLib::PoseData> output = HomLib::BarathVISAPP2016::get_haf(x, y, A, Fn);
            for (auto &p : output) poses->push_back(p);
            denormalize(poses, f0);
            return static_cast<int>(poses->size());
        }
        int minimal_sample_size() const {
            return 1;
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

class ThreePTSolver : public NoDistortionPointEstimator<ThreePTSolver> {
    public:
        static constexpr int minimal_sample_size_value = 3;
        static constexpr bool supports_overdetermined = true;
        Eigen::Matrix3d F_;
        explicit ThreePTSolver(const Eigen::Matrix3d &F) : F_(F) {}
        ThreePTSolver() : F_(Eigen::Matrix3d::Constant(std::numeric_limits<double>::quiet_NaN())) {}
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange & = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::BarathVISAPP2016::get_3pt(x, y, F_);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_,
                     std::vector<HomLib::PoseData> *poses) const {
            std::vector<Eigen::Vector2d> x = x_, y = y_;
            double f0;
            poses->clear();
            if (!normalize(x, y, &f0))
                return 0;
            Eigen::Matrix3d Fn = Eigen::Vector3d(f0, f0, 1).asDiagonal() * F_ * Eigen::Vector3d(f0, f0, 1).asDiagonal();
            std::vector<HomLib::PoseData> output = HomLib::BarathVISAPP2016::get_3pt(x, y, Fn);
            for (auto &p : output) poses->push_back(p);
            denormalize(poses, f0);
            return static_cast<int>(poses->size());
        }
    };

}
}

#endif  // INCLUDES_HOMLIB_GET_BARATH_VISAPP_2016_HPP_
