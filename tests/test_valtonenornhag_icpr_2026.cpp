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
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_valtonenornhag_icpr_2026.hpp"
#include "posedata.hpp"

TEST_CASE("Valtonen Ornhag ICPR 2026 - AFFINE") {
    
	std::vector<Eigen::Vector2d> p1 = {
		Eigen::Vector2d(-0.304372059145248, -0.174125682172663),
		Eigen::Vector2d(-0.308082306270758,  0.100342577802908)
	};

	std::vector<Eigen::Vector2d> p2 = {
		Eigen::Vector2d(0.596062339231747, 0.449851750693698),
		Eigen::Vector2d(0.284590553127076, 0.20125992440591)
	};

	std::vector<Eigen::Matrix2d> A = {
		(Eigen::Matrix2d() << -0.742737063351745, -1.25720051435581,
		                      1.29931162730836, -0.953395533218031).finished(),
		(Eigen::Matrix2d() << -0.763314941736403, -1.05311577674003,
		                      1.13911319830585, -0.839438880618645).finished()
	};

    std::vector<HomLib::PoseData> posedata = HomLib::ValtonenOrnhagICPR2026::get_affine(p1, p2, A, false);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 1);

    // Test distortion parameters
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(-0.15508855476035885).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.0).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << -0.389294613852058, -0.452105197122382, 0.0574163046746314,
 				 0.475012108793918, -0.355897552980124,  0.279279658894568,
 				-0.152511379317942,  0.121821752305461,  0.407920793089342;


    REQUIRE(posedata[0].homography.isApprox(expected, tol));

}


TEST_CASE("Valtonen Ornhag ICPR 2026 - ORI") {
    std::vector<Eigen::Vector2d> p1 = {
		Eigen::Vector2d(0.611733207361777, -0.100815813483642),
		Eigen::Vector2d(0.506500688761625, 0.22654451870654),
		Eigen::Vector2d(0.618492431042353, 0.490634389234033),
		Eigen::Vector2d(0.446207436531052, -0.221543454776215)
	};

	std::vector<Eigen::Vector2d> p2 = {
		Eigen::Vector2d(0.152378076324475, 0.471173894309057),
		Eigen::Vector2d(0.35249479911981, 0.604057206211103),
		Eigen::Vector2d(0.576429540272504, 0.546190419126011),
		Eigen::Vector2d(-0.0960825885072252, 0.57183541266364)
	};

	std::vector<Eigen::Vector2d> ori = {
		Eigen::Vector2d(2.82616873462069, 2.11084467936878),
		Eigen::Vector2d(2.76670870212448, 1.5791745485192),
		Eigen::Vector2d(3.03912230291483, 1.40242855410085),
		Eigen::Vector2d(2.54398570702913, 2.11604504644583)
	};

    std::vector<HomLib::PoseData> posedata = HomLib::ValtonenOrnhagICPR2026::get_ori(p1, p2, ori);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 1);

    // Test distortion parameters
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(-0.18221254202290402).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.0).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << 0.400486830270973,  0.357931239582734, -0.141724104171464,
 				0.114670003794941, 0.0864762777928815,  0.190330050981093,
 				0.795029931876192, 0.0101246929389779, 0.0222653484078536;

    REQUIRE(posedata[0].homography.isApprox(expected, tol));

}
