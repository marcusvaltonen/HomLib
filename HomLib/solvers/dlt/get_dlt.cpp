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

// Standard 4-point DLT (Direct Linear Transform) for homography estimation.
// R. Hartley and A. Zisserman, Multiple View Geometry in Computer Vision, Cambridge
// University Press, 2003.

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include "posedata.hpp"
#include "get_dlt.hpp"
#include "general_position.hpp"

namespace HomLib {
namespace DLT {

template <bool CheckSample>
std::vector<HomLib::PoseData> get(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool use_qr
) {
    std::vector<HomLib::PoseData> output;

    const int n = static_cast<int>(x.size());
    if (n < 4) {
        return output;
    }
    // A minimal sample: no solution if three points are collinear, or if the points would be on
    // both sides of the vanishing line of H (see general_position.hpp). Larger samples aren't
    // checked.
    if constexpr (CheckSample) {
        if (n == 4 && !HomLib::is_valid_sample(x.data(), y.data())) {
            return output;
        }
    }
    Eigen::Matrix<double, Eigen::Dynamic, 9> A(2 * n, 9);
    for (int i = 0; i < n; i++) {
        const double x0 = x[i](0), x1 = x[i](1);
        const double y0 = y[i](0), y1 = y[i](1);
        A.row(2 * i)     << 0, 0, 0, -x0, -x1, -1,  y1*x0,  y1*x1,  y1;
        A.row(2 * i + 1) << x0, x1, 1,  0,   0,  0, -y0*x0, -y0*x1, -y0;
    }

    Eigen::Matrix<double, 9, 1> h;
    if (use_qr && n == 4) {
        Eigen::HouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> qr(A.transpose());
        h = qr.householderQ() * Eigen::Matrix<double, 9, 1>::Unit(8);
    } else {
        Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, 9>> svd(A, Eigen::ComputeFullV);
        h = svd.matrixV().col(8);
    }

    Eigen::Matrix3d H = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(h.data());

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
                                                 const std::vector<Eigen::Vector2d> &, bool);
template std::vector<HomLib::PoseData> get<false>(const std::vector<Eigen::Vector2d> &,
                                                  const std::vector<Eigen::Vector2d> &, bool);

}  // namespace DLT
}  // namespace HomLib
