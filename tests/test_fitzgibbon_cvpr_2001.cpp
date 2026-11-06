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
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "posedata.hpp"

TEST_CASE("Fitzgibbon CVPR 2001") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.291910956401174,   0.827397315523355),
        Eigen::Vector2d( 0.998181103444196,   0.380888182504428),
        Eigen::Vector2d( 0.823771713189742,   0.834719771881603),
        Eigen::Vector2d(-0.352037947619727,  -0.327091812180983),
        Eigen::Vector2d(-0.645473006258758,  -0.378883379128205)
    };

    std::vector<Eigen::Vector2d> p2 = {
        Eigen::Vector2d(-0.380278821800649,  -0.340077827788294),
        Eigen::Vector2d( 0.528335392478284,  -0.680480654905349),
        Eigen::Vector2d( 0.957534928343294,   0.410935343442919),
        Eigen::Vector2d( 0.376937156586257,   0.882649958739618),
        Eigen::Vector2d( 0.462529355018110,   0.738720655679304)
    };

    std::vector<HomLib::PoseData> output = HomLib::FitzgibbonCVPR2001::get(p1, p2);
    
    REQUIRE(output.size() == 1);
    
    HomLib::PoseData posedata = output[0];

    double tol = 1e-12;

    // Test distortion parameter
    REQUIRE(posedata.distortion_parameter == Catch::Approx(-0.6554778901775142).margin(tol));

    // Test homography
    tol = 1e-7;
    Eigen::Matrix3d expected;
    expected << 0.07731617,  0.25714674, -0.13263548,
               -0.19494561,  0.39331832, -0.54704409,
                0.23663649, -0.11162851, -0.22216837;
    REQUIRE(posedata.homography.isApprox(expected, tol));
}

TEST_CASE("Fitzgibbon CVPR 2001 - eigensolver without convergence") {
    // A NaN coordinate stops the generalized eigensolver from converging: the solvers must return
    // no solution, without reading the eigenvalues and eigenvectors that Eigen didn't compute.
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(0.1, 0.2), Eigen::Vector2d(0.8, -0.3), Eigen::Vector2d(-0.5, 0.6),
        Eigen::Vector2d(0.4, 0.9), Eigen::Vector2d(-0.7, -0.4)
    };
    std::vector<Eigen::Vector2d> p2 = p1;
    p2[2](0) = std::numeric_limits<double>::quiet_NaN();

    REQUIRE(HomLib::FitzgibbonCVPR2001::get(p1, p2).empty());
    REQUIRE(HomLib::FitzgibbonCVPR2001::get_single(p1, p2).empty());
}
