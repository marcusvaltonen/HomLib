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
#include <algorithm>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_kukelova_cvpr_2015.hpp"
#include "posedata.hpp"

TEST_CASE("Kukelova CVPR 2015") {
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

    std::vector<HomLib::PoseData> posedata = HomLib::KukelovaCVPR2015::get(p1, p2, false);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 5);

    // Sort by distortion_parameter so the test is independent of eigenvalue ordering
    std::sort(posedata.begin(), posedata.end(),
              [](const HomLib::PoseData &a, const HomLib::PoseData &b) {
                  return a.distortion_parameter < b.distortion_parameter;
              });

    // Test distortion parameters (sorted ascending by distortion_parameter)
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(-2.136402668559706).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.6928831549898077).margin(tol));
    REQUIRE(posedata[1].distortion_parameter == Catch::Approx(-1.7814746086320201).margin(tol));
    REQUIRE(posedata[1].distortion_parameter2 == Catch::Approx(-2.079301697963529).margin(tol));
    REQUIRE(posedata[2].distortion_parameter == Catch::Approx(-1.1832137508476386).margin(tol));
    REQUIRE(posedata[2].distortion_parameter2 == Catch::Approx(-1.8809663034629707).margin(tol));
    REQUIRE(posedata[3].distortion_parameter == Catch::Approx(-0.6554778901775545).margin(tol));
    REQUIRE(posedata[3].distortion_parameter2 == Catch::Approx(-0.6554778901775485).margin(tol));
    REQUIRE(posedata[4].distortion_parameter == Catch::Approx(0.2537509800067458).margin(tol));
    REQUIRE(posedata[4].distortion_parameter2 == Catch::Approx(-1.9785613160596929).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << 0.34484818, 0.34901380, 0.35962541,
                0.57428832, 0.55155001, 0.50852640,
                0.99896668, 1.40954432, 1.14642718;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));

    expected << 0.16258481,  0.42614499,  0.33353189,
                0.34074170,  0.54233918,  0.50396602,
               -0.20919435, -0.52457951, -0.41537822;
    REQUIRE(posedata[1].homography.isApprox(expected, tol));

    expected << -0.07243639,  0.38036240,  0.27070323,
                 0.07137307,  0.35916405,  0.51141177,
                 0.04851674, -0.51620915, -0.47045628;
    REQUIRE(posedata[2].homography.isApprox(expected, tol));

    expected << -0.08070970, -0.26843331,  0.13845705,
                 0.20350208, -0.41058168,  0.57105471,
                -0.24702284,  0.11652806,  0.23191969;
    REQUIRE(posedata[3].homography.isApprox(expected, tol));

    expected << -1.84993859,  4.64801961,  0.72256438,
                -2.97548333,  4.41549851,  0.14438585,
                 2.57097121, -6.03637522, -0.79542495;
    REQUIRE(posedata[4].homography.isApprox(expected, tol));
}

TEST_CASE("Kukelova CVPR 2015 - range check of the estimator") {
    // A two-sided problem with k1 = k2 = -0.3 (the points of the Wadenbäck two-sided test)
    const std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d(-0.093972209335993,  0.0431909090094019),
        Eigen::Vector2d(-0.0969912533237537,  0.0171160938161422),
        Eigen::Vector2d(-0.123299663772998,  0.500859415499538),
        Eigen::Vector2d(-0.300610529524079,  0.477995136756032),
        Eigen::Vector2d( 0.143634675965547, 0.402159621931429)
    };
    const std::vector<Eigen::Vector2d> y = {
        Eigen::Vector2d(-0.0150485985127741,   0.504312576320796),
        Eigen::Vector2d(-0.0324017385815406,   0.525321116499528),
        Eigen::Vector2d(0.335015888737872, 0.182133325842435),
        Eigen::Vector2d(0.448124752741167, 0.313117448200805),
        Eigen::Vector2d(0.0822796392160207,  0.0746306773026531)
    };
    auto has_truth = [](const std::vector<HomLib::PoseData> &poses, double k) {
        for (const HomLib::PoseData &p : poses)
            if (std::abs(p.distortion_parameter - k) < 1e-5 * std::abs(k) && std::abs(p.distortion_parameter2 - k) < 1e-5 * std::abs(k))
                return true;
        return false;
    };
    auto norm2 = [](const std::vector<Eigen::Vector2d> &p) {   // f0^2 of the estimator
        double f0 = 0;
        for (const Eigen::Vector2d &q : p)
            f0 += q.norm();
        f0 /= p.size() * std::sqrt(2.0);
        return f0 * f0;
    };
    HomLib::KukelovaCVPR2015::SolverTwoSided solver;
    std::vector<HomLib::PoseData> all, poses;

    // Without a range
    solver.k_range = HomLib::DistortionRange();
    solver.estimate(x, y, &all);
    REQUIRE(has_truth(all, -0.3));

    // The default, [-1, 0.1] in the normalized coordinates (there, the true k is -0.3 f0^2)
    solver.k_range = HomLib::DistortionRange(-1.0, 0.1);
    solver.estimate(x, y, &poses);
    REQUIRE(has_truth(poses, -0.3));
    REQUIRE(poses.size() < all.size());
    for (const HomLib::PoseData &p : poses) {
        REQUIRE(p.distortion_parameter * norm2(x) >= -1.0 - 1e-12);
        REQUIRE(p.distortion_parameter * norm2(x) <= 0.1 + 1e-12);
        REQUIRE(p.distortion_parameter2 * norm2(x) >= -1.0 - 1e-12);
        REQUIRE(p.distortion_parameter2 * norm2(x) <= 0.1 + 1e-12);
    }

    // The same problem with the points scaled by s (k scaled by 1/s^2): the same solutions
    const double s = 1000.0;
    std::vector<Eigen::Vector2d> xs, ys;
    for (int i = 0; i < 5; i++) {
        xs.push_back(s * x[i]);
        ys.push_back(s * y[i]);
    }
    std::vector<HomLib::PoseData> scaled;
    solver.estimate(xs, ys, &scaled);
    REQUIRE(scaled.size() == poses.size());
    REQUIRE(has_truth(scaled, -0.3 / (s * s)));
    // A range that excludes the true k (in the normalized coordinates)
    solver.k_range = HomLib::DistortionRange(-0.2 * norm2(x), 0.1);
    solver.estimate(xs, ys, &scaled);
    REQUIRE_FALSE(has_truth(scaled, -0.3 / (s * s)));
}
