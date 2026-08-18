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
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "affine_pose_estimator.h"
#include "refinement.hpp"

namespace HomLib {
namespace BarathVISAPP2016 {
std::vector<HomLib::PoseData> get_affine(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Matrix2d> &A
);

class AffineSolver : public AffinePoseEstimator<AffineSolver> {
    public:
        AffineSolver() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::BarathVISAPP2016::get_affine(x, y, A);
            for (size_t i = 0; i < output.size(); i++) {
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

#endif  // INCLUDES_HOMLIB_GET_BARATH_VISAPP_2016_HPP_
