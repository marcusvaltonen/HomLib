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
#include <limits>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_nakano_icpr_2025.hpp"
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

TEST_CASE("Nakano ICPR 2025 - eigensolver without convergence") {
    // A NaN coordinate stops the eigensolver from converging: the solver must return no solution
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(0.1, 0.2), Eigen::Vector2d(0.8, -0.3), Eigen::Vector2d(-0.5, 0.6),
        Eigen::Vector2d(0.4, 0.9), Eigen::Vector2d(-0.7, -0.4)
    };
    std::vector<Eigen::Vector2d> p2 = p1;
    p2[2](0) = std::numeric_limits<double>::quiet_NaN();

    REQUIRE(HomLib::NakanoICPR2025::get(p1, p2, false).empty());
    REQUIRE(HomLib::NakanoICPR2025::get(p1, p2, true).empty());
}
