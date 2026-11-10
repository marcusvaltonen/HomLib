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

namespace HomLib {
namespace NakanoICPR2025 {
// N >= 5 points. For N == 5, uses QR for exact nullspace; for N > 5,
// uses eigendecomposition of M^T M for the approximate nullspace.
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check,
    const HomLib::DistortionRange &range = HomLib::DistortionRange()
);
class SolverSingleSided : public OneSidedLeftPointEstimator<SolverSingleSided> {
    public:
        static constexpr int minimal_sample_size_value = 5;
        bool extra_check = false;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::NakanoICPR2025::get(x, y, extra_check, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
    };
}
}

#endif  // INCLUDES_HOMLIB_GET_NAKANO_ICPR_2025_HPP_
