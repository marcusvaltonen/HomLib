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

// Adapted from TheiaSfM

#include <Eigen/Geometry>
#include <vector>
#include <limits>
#include <iostream>
#include "posedata.hpp"

namespace HomLib {
namespace KukelovaCVPR2015 {

using Vector6d = Eigen::Matrix<double, 6, 1>;
using Array51d = Eigen::Array<double, 5, 1>;
using Matrix68d = Eigen::Matrix<double, 6, 8>;
using Matrix62d = Eigen::Matrix<double, 6, 2>;
using Matrix65d = Eigen::Matrix<double, 6, 5>;

std::vector<HomLib::PoseData> get_6pt(
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    bool dist_equal
) {
    Matrix62d X;
    Matrix62d U;

    for (int i = 0; i < 6; ++i) {
        X.row(i) = y[i];
        U.row(i) = x[i];
    }

    Matrix68d M;
    Vector6d u2 = U.col(0).array().square() + U.col(1).array().square();

    M.col(0) = -X.col(1).array() * U.col(0).array();
    M.col(1) = -X.col(1).array() * U.col(1).array();
    M.col(2) = -X.col(1);
    M.col(3) = X.col(0).array() * U.col(0).array();
    M.col(4) = X.col(0).array() * U.col(1).array();
    M.col(5) = X.col(0);
    M.col(6) = -X.col(1).array() * u2.array();
    M.col(7) = X.col(0).array() * u2.array();

    Eigen::JacobiSVD<Matrix68d, Eigen::FullPivHouseholderQRPreconditioner> Svd1(
      M, Eigen::ComputeFullV);
    const Eigen::Matrix<double, 8, 8>& V1 = Svd1.matrixV();

    const double a = -V1(2, 6) * V1(7, 6) + V1(5, 6) * V1(6, 6);
    const double b = -V1(2, 6) * V1(7, 7) - V1(2, 7) * V1(7, 6) +
                   V1(5, 6) * V1(6, 7) + V1(5, 7) * V1(6, 6);
    const double c = -V1(2, 7) * V1(7, 7) + V1(5, 7) * V1(6, 7);
    const double d = b * b - 4.0 * a * c;

    int nsols = 0;
    Eigen::Vector2d rs;
    std::vector<HomLib::PoseData> output;

    if (std::abs(d) < 100.0 * std::numeric_limits<double>::epsilon()) {
        nsols = 1;
        rs(0) = (-b) / (2.0 * a);
    } else if (d > 0.0) {
        nsols = 2;
        double d2 = std::sqrt(d);
        rs(0) = (-b + d2) / (2.0 * a);
        rs(1) = (-b - d2) / (2.0 * a);
    } else {
        return output;
    }

    const Vector6d x2 = X.col(0).array().square() + X.col(1).array().square();
    Vector6d u3, r;
    Eigen::Matrix<double, 8, 1> n;
    Matrix65d T;
    T.col(0) = -M.col(3);
    T.col(1) = -M.col(4);

    for (int i = 0; i < nsols; i++) {
        n = rs(i) * V1.col(6) + V1.col(7);
        const double l2 = n(6) / n(2);

        u3 = u3.Ones() + l2 * u2;
        r = n(0) * U.col(0) + n(1) * U.col(1) + n(2) * u3;

        T.col(2) = -X.col(0).array() * u3.array();
        T.col(3) = x2.array() * r.array();
        T.col(4) = r;

        Eigen::JacobiSVD<Matrix65d> Svd2(T, Eigen::ComputeFullV);
        Eigen::Matrix<double, 5, 1> v2 = Svd2.matrixV().col(4);

        v2.head(4) /= v2(4);
        const double l1 = v2(3);

        //if (dist_equal && (std::abs(l1 - l2) > 0.05 * std::sqrt(l1 * l2))) {
        //    continue;
        //}

        // Package output
        HomLib::PoseData pd;
        pd.homography << n(0), n(1), n(2), n(3), n(4), n(5), v2(0), v2(1), v2(2); 
        if (dist_equal) {
            pd.distortion_parameter = 0.5 * (l1 + l2);
        } else {
            pd.distortion_parameter = l2;
            pd.distortion_parameter2 = l1;
        }
        output.push_back(pd);
    }

    return output;
}
}  // namespace KukelovaCVPR2015
}  // namespace HomLib
