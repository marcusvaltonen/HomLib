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

// 3PT method from Barath & Hajder, VISAPP 2016, Section 2.3.
// Homography from 3+ point correspondences and known Fundamental matrix.
// Substitutes [e2]x H = F into the DLT equations Hp1 ~ p2 to get a
// 2n x 3 system in (h31, h32, h33).
// Minimal: 3 points (overdetermined 6x3).

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "posedata.hpp"
#include "get_barath_visapp_2016.hpp"

namespace HomLib {
namespace BarathVISAPP2016 {

std::vector<HomLib::PoseData> get_3pt(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    const Eigen::Matrix3d &F
) {
    std::vector<HomLib::PoseData> output;
    if (x.size() < 3) return output;

    // Epipole e2 in image 2: F^T e2 = 0
    Eigen::JacobiSVD<Eigen::Matrix3d> svd_F(F.transpose(), Eigen::ComputeFullV);
    Eigen::Vector3d e2 = svd_F.matrixV().col(2);
    if (std::abs(e2(2)) < 1e-14) return output;
    e2 /= e2(2);
    const double ex = e2(0), ey = e2(1);

    const double f11 = F(0, 0), f12 = F(0, 1), f13 = F(0, 2);
    const double f21 = F(1, 0), f22 = F(1, 1), f23 = F(1, 2);

    const int n = static_cast<int>(x.size());
    Eigen::MatrixXd C(2 * n, 3);
    Eigen::VectorXd d(2 * n);

    for (int i = 0; i < n; i++) {
        const double u1 = x[i](0), v1 = x[i](1);
        const double u2 = y[i](0), v2 = y[i](1);
        const int row = 2 * i;

        // Substituting Eq. 4 into Hp1 ~ p2:
        // (ex - u2)(h31*u1 + h32*v1 + h33) = -(f21*u1 + f22*v1 + f23)
        C(row, 0) = (ex - u2) * u1;
        C(row, 1) = (ex - u2) * v1;
        C(row, 2) = (ex - u2);
        d(row) = -(f21 * u1 + f22 * v1 + f23);

        // (ey - v2)(h31*u1 + h32*v1 + h33) = (f11*u1 + f12*v1 + f13)
        C(row + 1, 0) = (ey - v2) * u1;
        C(row + 1, 1) = (ey - v2) * v1;
        C(row + 1, 2) = (ey - v2);
        d(row + 1) = f11 * u1 + f12 * v1 + f13;
    }

    Eigen::Vector3d h3 = C.colPivHouseholderQr().solve(d);

    // Recover full H from Eq. 4 (lambda = 1)
    Eigen::Matrix3d H;
    H(0, 0) = ex * h3(0) + f21;   H(0, 1) = ex * h3(1) + f22;   H(0, 2) = ex * h3(2) + f23;
    H(1, 0) = ey * h3(0) - f11;   H(1, 1) = ey * h3(1) - f12;   H(1, 2) = ey * h3(2) - f13;
    H(2, 0) = h3(0);              H(2, 1) = h3(1);              H(2, 2) = h3(2);

    if (H(2, 2) != 0.0) H /= H(2, 2);
    if (!H.allFinite()) return output;

    HomLib::PoseData pd;
    pd.homography = H;
    pd.distortion_parameter = 0.0;
    pd.distortion_parameter2 = 0.0;
    output.push_back(pd);
    return output;
}

}  // namespace BarathVISAPP2016
}  // namespace HomLib
