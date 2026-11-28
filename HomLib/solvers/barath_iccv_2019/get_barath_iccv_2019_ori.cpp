// Copyright (C) 2019 Czech Technical University.
// All rights reserved.
//
// Original 3-orientation (ORB) homography solver from Barath & Kukelova (ICCV 2019).
// Algorithm: 3n×9 linear system (2 point + 1 rotation per correspondence), SVD nullspace.
// Ported from the GC-RANSAC coefficient structure — no optimizations.

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "posedata.hpp"
#include "get_barath_iccv_2019.hpp"

namespace HomLib {
namespace BarathICCV2019 {

std::vector<HomLib::PoseData> get_ori(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    const std::vector<Eigen::Vector2d> &ori
) {
    std::vector<HomLib::PoseData> output;
    const int n = static_cast<int>(x.size());
    if (n < 3)
        return output;

    Eigen::MatrixXd coefficients(3 * n, 9);

    for (int i = 0; i < n; ++i) {
        const double u1 = x[i](0), v1 = x[i](1);
        const double u2 = y[i](0), v2 = y[i](1);
        const double s1 = std::sin(ori[i](0)), c1 = std::cos(ori[i](0));
        const double s2 = std::sin(ori[i](1)), c2 = std::cos(ori[i](1));

        const double c1s2 = c1 * s2;
        const double c1c2 = c1 * c2;
        const double s1s2 = s1 * s2;
        const double s1c2 = s1 * c2;

        coefficients.row(3 * i)     << 0, 0, 0, u1, v1, 1, -u1 * v2, -v1 * v2, -v2;
        coefficients.row(3 * i + 1) << u1, v1, 1, 0, 0, 0, -u1 * u2, -v1 * u2, -u2;
        coefficients.row(3 * i + 2) << -c1s2, -s1s2, 0, c1c2, s1c2, 0,
            u2 * c1s2 - v2 * c1c2, u2 * s1s2 - v2 * s1c2, 0;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(coefficients, Eigen::ComputeFullV);
    Eigen::VectorXd h = svd.matrixV().col(8);

    Eigen::Matrix3d H;
    H << h(0), h(1), h(2),
         h(3), h(4), h(5),
         h(6), h(7), h(8);

    if (H(2, 2) != 0.0)
        H /= H(2, 2);
    if (!H.allFinite())
        return output;

    HomLib::PoseData pd;
    pd.homography = H;
    pd.distortion_parameter = 0.0;
    pd.distortion_parameter2 = 0.0;
    output.push_back(pd);

    return output;
}

}  // namespace BarathICCV2019
}  // namespace HomLib
