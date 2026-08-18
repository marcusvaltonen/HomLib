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

#ifndef INCLUDES_HOMLIB_GET_NAKANO_ICPR_2025_HPP_
#define INCLUDES_HOMLIB_GET_NAKANO_ICPR_2025_HPP_

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "affine_pose_estimator.h"
#include "orientation_pose_estimator.h"
#include "pose_estimator.h"
#include "refinement.hpp"

namespace HomLib {
namespace NakanoICPR2025 {
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
);
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
std::vector<HomLib::PoseData> get_affine_no_dist(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Matrix2d> &A
);
class SolverSingleSided : public PoseEstimator<SolverSingleSided> {
    public:
        SolverSingleSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::NakanoICPR2025::get(x, y, extra_check);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, 0.0);  // One-sided
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter2);
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided(x, y, pose);
        }
    private:
        bool extra_check = false;
    };
class SolverSingleSidedRight : public PoseEstimator<SolverSingleSidedRight> {
    public:
        SolverSingleSidedRight() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::NakanoICPR2025::get(y, x, extra_check);  // Hack here
            for (size_t i = 0; i < output.size(); i++) {
                output[i].homography = output[i].homography.inverse();  //Hack here
                output[i].distortion_parameter = output[i].distortion_parameter2;  //Hack here
                output[i].distortion_parameter2 = 0.0;  //Hack here
                // std::cout << "H[" << i << "] = " << output[i].homography / output[i].homography(2,2) << std::endl;
                // std::cout << "k1[" << i << "] = " << output[i].distortion_parameter << std::endl;
                // std::cout << "k2[" << i << "] = " << output[i].distortion_parameter2 << std::endl;
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);  // One-sided
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = yu;
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided_right(x, y, pose);
        }
    private:
        bool extra_check = false;
    };

class AffineSolverSingleSided : public AffinePoseEstimator<AffineSolverSingleSided> {
    public:
        AffineSolverSingleSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::NakanoICPR2025::get_affine(x, y, A, extra_check);
            for (size_t i = 0; i < output.size(); i++) {
                // std::cout << "H[" << i << "] = " << output[i].homography / output[i].homography(2,2) << std::endl;
                // std::cout << "k1[" << i << "] = " << output[i].distortion_parameter << std::endl;
                // std::cout << "k2[" << i << "] = " << output[i].distortion_parameter2 << std::endl;
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 2;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);  // One-sided
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = yu;
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided_right(x, y, pose);
        }
    private:
        bool extra_check = false;
    };

class OrientationSolverSingleSided : public OrientationPoseEstimator<OrientationSolverSingleSided> {
    public:
        OrientationSolverSingleSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Vector2d> &ori, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::NakanoICPR2025::get_ori(x, y, ori);
            for (size_t i = 0; i < output.size(); i++) {
                // std::cout << "H[" << i << "] = " << output[i].homography / output[i].homography(2,2) << std::endl;
                // std::cout << "k1[" << i << "] = " << output[i].distortion_parameter << std::endl;
                // std::cout << "k2[" << i << "] = " << output[i].distortion_parameter2 << std::endl;
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 4;  // Degenerates ?
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);  // One-sided
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = yu;
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided_right(x, y, pose);
        }
    };

class AffineSolverNoDist : public AffinePoseEstimator<AffineSolverNoDist> {
    public:
        AffineSolverNoDist() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::NakanoICPR2025::get_affine_no_dist(x, y, A);
            for (size_t i = 0; i < output.size(); i++) {
                // std::cout << "H[" << i << "] = " << output[i].homography / output[i].homography(2,2) << std::endl;
                // std::cout << "k1[" << i << "] = " << output[i].distortion_parameter << std::endl;
                // std::cout << "k2[" << i << "] = " << output[i].distortion_parameter2 << std::endl;
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 2;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = xd;  // No dist
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = yu;  // No dist
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_no_dist(x, y, pose);
        }
    };

}
}

#endif  // INCLUDES_HOMLIB_GET_NAKANO_ICPR_2025_HPP_
