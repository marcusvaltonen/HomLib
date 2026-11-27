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

// HAF method from Barath & Hajder, VISAPP 2016, Section 2.3.
// Homography from Affine transformation and known Fundamental matrix.
// Substitutes [e2]x H = F into the affine equations to get an
// inhomogeneous 4n x 3 system in (h31, h32, h33).
// Minimal: 1 AC (overdetermined 4x3).

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "posedata.hpp"
#include "get_barath_visapp_2016.hpp"

namespace HomLib {
namespace BarathVISAPP2016 {

std::vector<HomLib::PoseData> get_haf(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    const std::vector<Eigen::Matrix2d> &Aff,
    const Eigen::Matrix3d &F
) {
    std::vector<HomLib::PoseData> output;
    if (x.empty()) return output;

    // Epipole e2 in image 2: F^T e2 = 0
    Eigen::JacobiSVD<Eigen::Matrix3d> svd_F(F.transpose(), Eigen::ComputeFullV);
    Eigen::Vector3d e2 = svd_F.matrixV().col(2);
    if (std::abs(e2(2)) < 1e-14) return output;
    e2 /= e2(2);
    const double ex = e2(0), ey = e2(1);

    const double f11 = F(0, 0), f12 = F(0, 1), f13 = F(0, 2);
    const double f21 = F(1, 0), f22 = F(1, 1), f23 = F(1, 2);

    const int n = static_cast<int>(x.size());
    Eigen::MatrixXd C(4 * n, 3);
    Eigen::VectorXd d(4 * n);

    for (int i = 0; i < n; i++) {
        const double u1 = x[i](0), v1 = x[i](1);
        const double u2 = y[i](0), v2 = y[i](1);
        const double a11 = Aff[i](0, 0), a12 = Aff[i](0, 1);
        const double a21 = Aff[i](1, 0), a22 = Aff[i](1, 1);
        const int row = 4 * i;

        // Eq. 8 from the paper (lambda = 1)
        C(row, 0) = a11 * u1 + u2 - ex;  C(row, 1) = a11 * v1;                C(row, 2) = a11;
        d(row) = f21;

        C(row + 1, 0) = a12 * u1;         C(row + 1, 1) = a12 * v1 + u2 - ex;  C(row + 1, 2) = a12;
        d(row + 1) = f22;

        C(row + 2, 0) = a21 * u1 + v2 - ey;  C(row + 2, 1) = a21 * v1;         C(row + 2, 2) = a21;
        d(row + 2) = -f11;

        C(row + 3, 0) = a22 * u1;         C(row + 3, 1) = a22 * v1 + v2 - ey;  C(row + 3, 2) = a22;
        d(row + 3) = -f12;
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
