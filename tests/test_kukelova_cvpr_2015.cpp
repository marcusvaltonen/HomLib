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
#include <vector>
#include <catch2/catch.hpp>
#include "get_kukelova_cvpr_2015.hpp"
#include "posedata.hpp"

TEST_CASE("Kukelova CVPR 2015") {
    Eigen::MatrixXd p1(2, 5);
    Eigen::MatrixXd p2(2, 5);
    p1 << -0.291910956401174,   0.998181103444196,   0.823771713189742,  -0.352037947619727,  -0.645473006258758,
           0.827397315523355,   0.380888182504428,   0.834719771881603,  -0.327091812180983,  -0.378883379128205;

    p2 << -0.380278821800649,   0.528335392478284,   0.957534928343294,   0.376937156586257,   0.462529355018110,
          -0.340077827788294,  -0.680480654905349,   0.410935343442919,   0.882649958739618,   0.738720655679304;

    std::vector<HomLib::PoseData> posedata = HomLib::KukelovaCVPR2015::get(p1, p2);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 5);

    // Test distortion parameters
    REQUIRE(posedata[0].distortion_parameter == Approx(0.2537509800067458).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Approx(-1.9785613160596929).margin(tol));
    REQUIRE(posedata[1].distortion_parameter == Approx(-1.1832137508476386).margin(tol));
    REQUIRE(posedata[1].distortion_parameter2 == Approx(-1.8809663034629707).margin(tol));
    REQUIRE(posedata[2].distortion_parameter == Approx(-1.7814746086320201).margin(tol));
    REQUIRE(posedata[2].distortion_parameter2 == Approx(-2.079301697963529).margin(tol));
    REQUIRE(posedata[3].distortion_parameter == Approx(-2.136402668559706).margin(tol));
    REQUIRE(posedata[3].distortion_parameter2 == Approx(0.6928831549898077).margin(tol));
    REQUIRE(posedata[4].distortion_parameter == Approx(-0.6554778901775545).margin(tol));
    REQUIRE(posedata[4].distortion_parameter2 == Approx(-0.6554778901775485).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << -1.84993859,  4.64801961,  0.72256438,
                -2.97548333,  4.41549851,  0.14438585,
                 2.57097121, -6.03637522, -0.79542495;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));

    expected << -0.07243639,  0.38036240,  0.27070323,
                 0.07137307,  0.35916405,  0.51141177,
                 0.04851674, -0.51620915, -0.47045628;
    REQUIRE(posedata[1].homography.isApprox(expected, tol));

    expected << 0.16258481,  0.42614499,  0.33353189,
                0.34074170,  0.54233918,  0.50396602,
               -0.20919435, -0.52457951, -0.41537822;
    REQUIRE(posedata[2].homography.isApprox(expected, tol));

    expected << 0.34484818, 0.34901380, 0.35962541,
                0.57428832, 0.55155001, 0.50852640,
                0.99896668, 1.40954432, 1.14642718;
    REQUIRE(posedata[3].homography.isApprox(expected, tol));

    expected << -0.08070970, -0.26843331,  0.13845705,
                 0.20350208, -0.41058168,  0.57105471,
                -0.24702284,  0.11652806,  0.23191969;
    REQUIRE(posedata[4].homography.isApprox(expected, tol));
}
