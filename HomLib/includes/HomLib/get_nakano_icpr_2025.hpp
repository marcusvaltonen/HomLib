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
#include "pose_estimator.h"
#include "refinement.hpp"
#include <iostream>
namespace HomLib {
namespace NakanoICPR2025 {
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
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
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter);
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided(x, y, pose);
        }
    private:
        bool extra_check = false;
    };
}
}

#endif  // INCLUDES_HOMLIB_GET_NAKANO_ICPR_2025_HPP_
