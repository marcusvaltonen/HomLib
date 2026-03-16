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

#ifndef INCLUDES_HOMLIB_GET_KUKELOVA_CVPR_2015_HPP_
#define INCLUDES_HOMLIB_GET_KUKELOVA_CVPR_2015_HPP_

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "pose_estimator.h"
#include "refinement.hpp"

namespace HomLib {
namespace KukelovaCVPR2015 {
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool dist_equal
);
std::vector<HomLib::PoseData> get_6pt(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool dist_equal
);
class SolverTwoSidedEqual : public PoseEstimator<SolverTwoSidedEqual> {
    public:
        SolverTwoSidedEqual() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get(x, y, dist_equal);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter);
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_twosided_equal(x, y, pose);
        }
    private:
        bool dist_equal = true;
    };
class SolverTwoSidedEqual6Pt : public PoseEstimator<SolverTwoSidedEqual6Pt> {
    public:
        SolverTwoSidedEqual6Pt() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get_6pt(x, y, dist_equal);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 6;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter);
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_twosided_equal(x, y, pose);
        }
    private:
        bool dist_equal = true;
    };
class SolverTwoSided : public PoseEstimator<SolverTwoSided> {
    public:
        SolverTwoSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get(x, y, dist_equal);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter2);
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_twosided(x, y, pose);
        }
    private:
        bool dist_equal = false;
    };
class SolverTwoSided6Pt : public PoseEstimator<SolverTwoSided6Pt> {
    public:
        SolverTwoSided6Pt() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get_6pt(x, y, dist_equal);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 6;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter2);
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_twosided(x, y, pose);
        }
    private:
        bool dist_equal = false;
    };
}
}

#endif  // INCLUDES_HOMLIB_GET_KUKELOVA_CVPR_2015_HPP_
