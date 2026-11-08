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

// S. Cai, Z. Wu, L. Guo, J. Wang, S. Zhang, J. Yan, S. Shen, Fast and interpretable 2d homography decomposition:
// Similarity-kernel-similarity and affine-core-affine transformations. PAMI 2025.
//
// Reference implementation: PoseLib (V. Larsson et al.)
// https://github.com/PoseLib/PoseLib/blob/master/PoseLib/solvers/homography_4pt.cc

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include "posedata.hpp"
#include "get_cai_pami_2025.hpp"
#include "general_position.hpp"

namespace HomLib {
namespace CaiPAMI2025 {

template <bool CheckSample>
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y
) {
    std::vector<HomLib::PoseData> output;
    if (x.size() < 4) {
        return output;
    }
    // No solution if three points are collinear, or if the points would be on both sides of the
    // vanishing line of H (see general_position.hpp)
    if constexpr (CheckSample) {
        if (!HomLib::is_valid_sample(x.data(), y.data())) {
            return output;
        }
    }

    // Source plane: differences from point 0
    const double M1N1_x = x[1](0) - x[0](0), M1P1_x = x[2](0) - x[0](0), M1Q1_x = x[3](0) - x[0](0);
    const double M1N1_y = x[1](1) - x[0](1), M1P1_y = x[2](1) - x[0](1), M1Q1_y = x[3](1) - x[0](1);
    const double fA1 = M1N1_x * M1P1_y - M1N1_y * M1P1_x;
    const double Q3_x = M1P1_y * M1Q1_x - M1P1_x * M1Q1_y;
    const double Q3_y = M1N1_x * M1Q1_y - M1N1_y * M1Q1_x;

    // Target plane: differences from point 0
    const double M2N2_x = y[1](0) - y[0](0), M2P2_x = y[2](0) - y[0](0), M2Q2_x = y[3](0) - y[0](0);
    const double M2N2_y = y[1](1) - y[0](1), M2P2_y = y[2](1) - y[0](1), M2Q2_y = y[3](1) - y[0](1);
    const double fA2 = M2N2_x * M2P2_y - M2N2_y * M2P2_x;
    const double Q4_x = M2P2_y * M2Q2_x - M2P2_x * M2Q2_y;
    const double Q4_y = M2N2_x * M2Q2_y - M2N2_y * M2Q2_x;

    // Core transformation H_C
    const double tt1 = fA1 - Q3_x - Q3_y;
    const double C11 = Q3_y * Q4_x * tt1;
    const double C22 = Q3_x * Q4_y * tt1;
    const double C33 = Q3_x * Q3_y * (fA2 - Q4_x - Q4_y);
    const double C31 = C11 - C33;
    const double C32 = C22 - C33;

    // H1 = H_A2_inv * H_C (upper-left 2x2)
    const double tt3 = y[0](0) * C33;
    const double tt4 = y[0](1) * C33;
    const double H1_11 = y[1](0) * C11 - tt3;
    const double H1_12 = y[2](0) * C22 - tt3;
    const double H1_21 = y[1](1) * C11 - tt4;
    const double H1_22 = y[2](1) * C22 - tt4;

    // H = H1 * H_A1
    Eigen::Matrix<double, 9, 1> h;
    h[0] = H1_11 * M1P1_y - H1_12 * M1N1_y;
    h[1] = H1_12 * M1N1_x - H1_11 * M1P1_x;
    h[3] = H1_21 * M1P1_y - H1_22 * M1N1_y;
    h[4] = H1_22 * M1N1_x - H1_21 * M1P1_x;
    h[6] = C31 * M1P1_y - C32 * M1N1_y;
    h[7] = C32 * M1N1_x - C31 * M1P1_x;
    h[2] = tt3 * fA1 - h[0] * x[0](0) - h[1] * x[0](1);
    h[5] = tt4 * fA1 - h[3] * x[0](0) - h[4] * x[0](1);
    h[8] = C33 * fA1 - h[6] * x[0](0) - h[7] * x[0](1);

    Eigen::Matrix3d H = Eigen::Map<const Eigen::Matrix3d>(h.data()).transpose();

    H.normalize();
    if (!H.allFinite() || std::abs(H.determinant()) < 1e-8) {
        return output;
    }

    HomLib::PoseData pd;
    pd.homography = H;
    output.push_back(pd);

    return output;
}

template std::vector<HomLib::PoseData> get<true>(const std::vector<Eigen::Vector2d> &,
                                                 const std::vector<Eigen::Vector2d> &);
template std::vector<HomLib::PoseData> get<false>(const std::vector<Eigen::Vector2d> &,
                                                  const std::vector<Eigen::Vector2d> &);

}  // namespace CaiPAMI2025
}  // namespace HomLib
