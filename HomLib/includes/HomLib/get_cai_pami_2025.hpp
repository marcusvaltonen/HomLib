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

#ifndef INCLUDES_HOMLIB_GET_CAI_PAMI_2025_HPP_
#define INCLUDES_HOMLIB_GET_CAI_PAMI_2025_HPP_

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "pose_estimator.h"
#include "distortion_range.h"
#include "refinement.hpp"

namespace HomLib {
namespace CaiPAMI2025 {

// 4-point homography using the ACA (Affine-Core-Affine) decomposition.
// No distortion estimation; returns 0 or 1 solution.
// CheckSample = false skips the check of the sample (is_valid_sample): faster, but an invalid
// sample can then give a wrong homography.
template <bool CheckSample = true>
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y
);

class Solver : public PoseEstimator<Solver> {
    public:
        static constexpr int minimal_sample_size_value = 4;
        static constexpr DistortionCase DISTORTION_CASE = NO_DISTORTION;
        Solver() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange & = HomLib::DistortionRange()) const {   // no distortion (k = 0 passes)
            std::vector<HomLib::PoseData> output = HomLib::CaiPAMI2025::get(x, y);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 4;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData, const Eigen::Vector2d &xd) const {
            return xd;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData, const Eigen::Vector2d &yu) const {
            return yu;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_no_dist(x, y, pose);
        }
    };

}  // namespace CaiPAMI2025
}  // namespace HomLib

#endif  // INCLUDES_HOMLIB_GET_CAI_PAMI_2025_HPP_
