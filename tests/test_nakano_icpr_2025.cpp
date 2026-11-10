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

#include <Eigen/Dense>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_nakano_icpr_2025.hpp"
#include "radial.hpp"
#include "posedata.hpp"

TEST_CASE("Nakano ICPR 2025") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.514683699631253, -0.0622906242414419),
        Eigen::Vector2d(-0.216333891824487,  0.236104964783804),
        Eigen::Vector2d( 0.0149270190286551, -0.466167744934018),
        Eigen::Vector2d(-0.299909368807735, -0.631181204140851),
        Eigen::Vector2d(-0.0667164585547235, 0.113775620793983)
    };

    std::vector<Eigen::Vector2d> p2 = {
        Eigen::Vector2d(-0.173489966290854,  0.161874840375985),
        Eigen::Vector2d(-0.276098980665827,  0.528540478431576),
        Eigen::Vector2d( 0.389704168705386,  0.291358438524245),
        Eigen::Vector2d( 0.307091630022996,  0.0216079362297207),
        Eigen::Vector2d(-0.0999459278021077, 0.583861165740761)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::NakanoICPR2025::get(p1, p2, false);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 1);

    // Test distortion parameters
    // One-sided (left): solver undistorts y, so distortion is in distortion_parameter2.
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(0.0).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(-0.3).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << -0.966294723999024,  1.424868744982000, -0.0570986060219543,
                -1.118034960806880, -0.916112342444385, -0.9604495174992110,
                 0.851682151867801,  0.527896230136481, -1.5204980422764300;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));

}

// ============================================================================
// Overdetermined: N > 5
// ============================================================================

static void generate_one_sided_left(
    const Eigen::Matrix3d& H, double k2, int n,
    std::vector<Eigen::Vector2d>* x, std::vector<Eigen::Vector2d>* y,
    unsigned seed = 42)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> coord(-0.8, 0.8);
    x->resize(n);
    y->resize(n);
    for (int i = 0; i < n; i++) {
        (*x)[i] = Eigen::Vector2d(coord(rng), coord(rng));
        Eigen::Vector3d tmp = H * (*x)[i].homogeneous();
        Eigen::Vector2d yu = tmp.hnormalized();
        Eigen::Vector2d yd;
        REQUIRE(HomLib::try_radialdistort(yu, k2, yd));
        (*y)[i] = yd;
    }
}

TEST_CASE("Nakano - overdetermined N > 5") {
    Eigen::Matrix3d H;
    H << 0.9, 0.1, 0.05,
        -0.1, 0.85, 0.03,
         0.02, -0.01, 1.0;
    double k2 = -0.15;

    std::vector<Eigen::Vector2d> x, y;
    generate_one_sided_left(H, k2, 20, &x, &y);

    auto output = HomLib::NakanoICPR2025::get(x, y, true);
    REQUIRE(output.size() == 1);
    REQUIRE(output[0].distortion_parameter2 == Catch::Approx(k2).margin(1e-3));

    Eigen::Matrix3d H_est = output[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_gt = H / H(2, 2);
    double max_reproj = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        Eigen::Vector2d yu_est = (H_est * x[i].homogeneous()).hnormalized();
        Eigen::Vector2d yu_gt = (H_gt * x[i].homogeneous()).hnormalized();
        max_reproj = std::max(max_reproj, (yu_est - yu_gt).norm());
    }
    REQUIRE(max_reproj < 1e-3);
}

TEST_CASE("Nakano - fewer than 5 returns empty") {
    std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d(0.1, 0.2), Eigen::Vector2d(0.3, 0.4),
        Eigen::Vector2d(0.5, 0.6), Eigen::Vector2d(0.7, 0.8)
    };
    REQUIRE(HomLib::NakanoICPR2025::get(x, x, false).empty());
}

TEST_CASE("Nakano - minimal N=5 regression") {
    Eigen::Matrix3d H;
    H << 0.9, 0.1, 0.05,
        -0.1, 0.85, 0.03,
         0.02, -0.01, 1.0;
    double k2 = -0.15;

    std::vector<Eigen::Vector2d> x, y;
    generate_one_sided_left(H, k2, 5, &x, &y);

    auto output = HomLib::NakanoICPR2025::get(x, y, true);
    REQUIRE(output.size() >= 1);
    REQUIRE(output[0].distortion_parameter2 == Catch::Approx(k2).margin(1e-6));
}
