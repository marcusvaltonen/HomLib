// Copyright (c) 2026 Frédéric Devernay
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

// Tests of the general-position checks (general_position.hpp) and of the solvers that use them.

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "general_position.hpp"
#include "get_cai_pami_2025.hpp"
#include "get_dlt.hpp"
#include "get_guo_prl_2019.hpp"
#include "get_wadenback_3dv_2026.hpp"

namespace {

Eigen::Vector2d apply(const Eigen::Matrix3d &H, const Eigen::Vector2d &x) {
    return (H * x.homogeneous()).hnormalized();
}

const Eigen::Matrix3d kH =
    (Eigen::Matrix3d() << 1.1, 0.2, 0.05, -0.1, 0.9, -0.03, 0.2, -0.1, 1.0).finished();

// A point of the division model, distorted with k, whose undistorted point is u
Eigen::Vector2d distort(const Eigen::Vector2d &u, double k) {
    return u * (2 / (1 + std::sqrt(1 - 4 * k * u.squaredNorm())));
}

}  // namespace

TEST_CASE("General position - four points") {
    Eigen::Vector2d x[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    Eigen::Vector2d y[4];
    for (int i = 0; i < 4; i++)
        y[i] = apply(kH, x[i]);
    REQUIRE(HomLib::is_valid_sample(x, y));

    // A projective homography with det H < 0 flips every triangle
    const Eigen::Matrix3d F =
        (Eigen::Matrix3d() << -1.0, 0.1, 0.2, 0.05, 1.0, 0.0, 0.3, 0.1, 1.0).finished();
    REQUIRE(F.determinant() < 0);
    for (int i = 0; i < 4; i++)
        y[i] = apply(F, x[i]);
    REQUIRE(HomLib::is_valid_sample(x, y));

    // H maps these points exactly, but its vanishing line (u = -1/2) separates them: some points
    // would be behind a camera
    const Eigen::Matrix3d V = (Eigen::Matrix3d() << 1, 0, 0, 0, 1, 0, 1, 0, 0.5).finished();
    Eigen::Vector2d u[4] = {{0, 0}, {1, 0}, {-1, 1}, {0, 1}};
    for (int i = 0; i < 4; i++)
        y[i] = apply(V, u[i]);
    REQUIRE_FALSE(HomLib::is_valid_sample(u, y));

    // Three collinear points, two equal ones, or a NaN
    for (int i = 0; i < 4; i++)
        y[i] = apply(kH, x[i]);
    Eigen::Vector2d c[4] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}};
    REQUIRE_FALSE(HomLib::is_valid_sample(c, y));
    REQUIRE_FALSE(HomLib::is_valid_sample(y, c));
    Eigen::Vector2d d[4] = {{0, 0}, {1, 0}, {1, 0}, {0, 1}};
    REQUIRE_FALSE(HomLib::is_valid_sample(d, y));
    Eigen::Vector2d e[4] = {{0, 0}, {1, 0}, {NAN, 1}, {0, 1}};
    REQUIRE_FALSE(HomLib::is_valid_sample(e, y));
}

TEST_CASE("General position - choosing four of five points") {
    // Points 0, 1, 2 are collinear: the choice leaves one of them out
    Eigen::Vector2d x[5] = {{0, 0}, {0.5, 0.5}, {1, 1}, {1, 0}, {0, 1}};
    Eigen::Vector2d y[5];
    for (int i = 0; i < 5; i++)
        y[i] = apply(kH, x[i]);
    REQUIRE(HomLib::select_valid_sample(x, y));
    REQUIRE(HomLib::is_valid_sample(x, y));
    for (int i = 0; i < 5; i++)
        REQUIRE((apply(kH, x[i]) - y[i]).norm() < 1e-12);   // still matched after the reordering

    // Every choice of 4 points has 3 collinear points
    Eigen::Vector2d c[5] = {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {0, 1}};
    for (int i = 0; i < 5; i++)
        y[i] = apply(kH, c[i]);
    REQUIRE_FALSE(HomLib::select_valid_sample(c, y));
}

TEST_CASE("General position - the 4-point solvers return nothing for three collinear points") {
    // Points 1, 2, 3 are collinear (with points 0, 1, 2, Guo's formula gives H = 0, which the
    // determinant test already rejected)
    const std::vector<Eigen::Vector2d> x = {{1, 0}, {0, 0}, {0.5, 0.5}, {1, 1}};
    std::vector<Eigen::Vector2d> y;
    for (const Eigen::Vector2d &p : x)
        y.push_back(apply(kH, p));
    REQUIRE(HomLib::GuoPRL2019::get(x, y).empty());       // was a NaN homography
    REQUIRE(HomLib::CaiPAMI2025::get(x, y).empty());
    REQUIRE(HomLib::DLT::get(x, y, false).empty());       // was an arbitrary homography
}

TEST_CASE("General position - Wadenback 3DV 2026 with three collinear undistorted points") {
    // Noise-free, k = -0.2 in both images; the undistorted points 0, 1, 2 are collinear, so the
    // first four points don't determine H
    const double k = -0.2;
    const std::vector<Eigen::Vector2d> u = {
        {-0.4, -0.2}, {0.0, 0.1}, {0.4, 0.4}, {0.5, -0.4}, {-0.3, 0.5}};
    std::vector<Eigen::Vector2d> x, y;
    for (const Eigen::Vector2d &p : u) {
        x.push_back(distort(p, k));
        y.push_back(distort(apply(kH, p), k));
    }
    for (bool extra_check : {false, true}) {
        std::vector<HomLib::PoseData> sols =
            HomLib::Wadenback3DV2026::get_double_sided_equal(x, y, extra_check);
        bool found = false;
        for (const HomLib::PoseData &s : sols) {
            if (std::abs(s.distortion_parameter - k) > 1e-9)
                continue;
            double err = 0;
            for (const Eigen::Vector2d &p : u)
                err = std::max(err, (apply(s.homography, p) - apply(kH, p)).norm());
            found = found || err < 1e-9;
        }
        CAPTURE(extra_check);
        REQUIRE(found);
    }
}
