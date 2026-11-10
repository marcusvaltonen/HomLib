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
#include "get_fitzgibbon_cvpr_2001.hpp"

namespace HomLib {
namespace KukelovaCVPR2015 {
// With a range: keeping only the solutions with distortion parameters in the range (checked
// before computing their homography)
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool dist_equal,
    const HomLib::DistortionRange &range = HomLib::DistortionRange()
);
// With a range: keeping only the solutions with distortion parameters in the range (checked
// before computing their homography)
std::vector<HomLib::PoseData> get_6pt(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool dist_equal,
    const HomLib::DistortionRange &range = HomLib::DistortionRange()
);
class SolverTwoSidedEqual : public TwoSidedEqualPointEstimator<SolverTwoSidedEqual> {
    public:
        static constexpr int minimal_sample_size_value = 5;
        SolverTwoSidedEqual() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get(x, y, dist_equal, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
    private:
        bool dist_equal = true;
    };
class SolverTwoSidedEqual6Pt : public TwoSidedEqualPointEstimator<SolverTwoSidedEqual6Pt> {
    public:
        static constexpr int minimal_sample_size_value = 6;
        SolverTwoSidedEqual6Pt() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get_6pt(x, y, dist_equal, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
    private:
        bool dist_equal = true;
    };
class SolverTwoSided : public TwoSidedPointEstimator<SolverTwoSided> {
    public:
        static constexpr int minimal_sample_size_value = 5;
        SolverTwoSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get(x, y, dist_equal, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
    private:
        bool dist_equal = false;
    };
class SolverTwoSided6Pt : public TwoSidedPointEstimator<SolverTwoSided6Pt> {
    public:
        static constexpr int minimal_sample_size_value = 6;
        SolverTwoSided6Pt() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::KukelovaCVPR2015::get_6pt(x, y, dist_equal, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
    private:
        bool dist_equal = false;
    };
}
}

#endif  // INCLUDES_HOMLIB_GET_KUKELOVA_CVPR_2015_HPP_
