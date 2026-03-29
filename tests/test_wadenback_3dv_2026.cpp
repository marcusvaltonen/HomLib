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
#include "get_wadenback_3dv_2026.hpp"
#include "posedata.hpp"


TEST_CASE("Wadenback 3DV 2026 - one-sided") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d( 0.242104535827652, -0.296297354299645),
        Eigen::Vector2d( 0.319776775526576, -0.114580331777873),
        Eigen::Vector2d(-0.122936754370192, -0.532655969828257),
        Eigen::Vector2d( 0.303141152380151, -0.295669384794041),
        Eigen::Vector2d( 0.446905125282507, -0.528127859697035)
    };

    std::vector<Eigen::Vector2d> p2 = {
        Eigen::Vector2d(-0.562787932824692, -0.0769625591786919),
        Eigen::Vector2d(-0.612813481767545, -0.0311462459576149),
        Eigen::Vector2d(-0.577247052139350, -0.252034297773712),
        Eigen::Vector2d(-0.544031097770073, -0.0440402280844913),
        Eigen::Vector2d(-0.424602042497788,  0.0211895908990797)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::Wadenback3DV2026::get_one_sided(p1, p2, false);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 3);

    // Test distortion parameters
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(5.11091).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.0).margin(tol));
    REQUIRE(posedata[1].distortion_parameter == Catch::Approx(2.73767).margin(tol));
    REQUIRE(posedata[1].distortion_parameter2 == Catch::Approx(0.0).margin(tol));
    REQUIRE(posedata[2].distortion_parameter == Catch::Approx(-0.3).margin(tol));
    REQUIRE(posedata[2].distortion_parameter2 == Catch::Approx(0.0).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected <<  3.24632054513494e-05, -4.95598175719004e-05, -1.87243913120627e-05,
                -8.71803482662030e-07, -3.24711594819769e-06, -2.28713196044006e-07,
                -0.000146151989535363,  0.000231387475606354,  8.59648988309500e-05;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));
    
    expected <<  7.88700481939429e-05, -0.000115003527994859, -5.31700374887011e-05,
                 6.35935433111556e-06, -6.84986213748709e-06, -3.56922455715305e-06,
                -0.000263190991107698,  0.000381884387848442,  0.000176871066503982;
    REQUIRE(posedata[1].homography.isApprox(expected, tol));
    
    expected <<  0.000227753625098076, -7.00643123791218e-05, -0.000525931983848800,
                 0.000429490139445578,  3.06112333501461e-05, -0.000156454396552730,
                 0.000133561563704977, -0.000486333541936380,  0.000545809522795618;
    REQUIRE(posedata[2].homography.isApprox(expected, tol));

}



/*
Wadenback 2025 (4.5 pt)
Not found: homography =
-13.6284238835806,  14.1792115247234, -2.18751115985823;
-13.5032381875633, -12.3737062846023,  8.92508296873142;
 4.57310999703002,  6.95986054985586,  17.8251082781653
distortion_parameter = -0.3
distortion_parameter2 = -0.3
x1 =
-0.093972209335993;
0.0431909090094019
x2 =
-0.0150485985127741;
  0.504312576320796
x1 =
-0.0969912533237537;
 0.0171160938161422
x2 =
-0.0324017385815406;
  0.525321116499528
x1 =
-0.123299663772998;
 0.500859415499538
x2 =
0.335015888737872;
0.182133325842435
x1 =
-0.300610529524079;
 0.477995136756032
x2 =
0.448124752741167;
0.313117448200805
x1 =
0.143634675965547;
0.402159621931429
x2 =
0.0822796392160207;
0.0746306773026531
pd[0].homography =
-5.37070541971519e-06,  5.58776046545164e-06, -8.62056987828971e-07;
-5.32137209241095e-06, -4.87624482275744e-06,  3.51720726336518e-06;
 1.80217660206107e-06,  2.74275008576624e-06,  7.02453977468668e-06
pd[0].distortion_parameter = -0.3
pd[0].distortion_parameter2 = -0.3

*/




TEST_CASE("Wadenback 3DV 2026 - two-sided equal") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.093972209335993,  0.0431909090094019),
        Eigen::Vector2d(-0.0969912533237537,  0.0171160938161422),
        Eigen::Vector2d(-0.123299663772998,  0.500859415499538),
        Eigen::Vector2d(-0.300610529524079,  0.477995136756032),
        Eigen::Vector2d( 0.143634675965547, 0.402159621931429)
    };

    std::vector<Eigen::Vector2d> p2 = {
        Eigen::Vector2d(-0.0150485985127741,   0.504312576320796),
        Eigen::Vector2d(-0.0324017385815406,   0.525321116499528),
        Eigen::Vector2d(0.335015888737872, 0.182133325842435),
        Eigen::Vector2d(0.448124752741167, 0.313117448200805),
        Eigen::Vector2d(0.0822796392160207,  0.0746306773026531)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::Wadenback3DV2026::get_double_sided_equal(p1, p2, true);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 1);

    // Test distortion parameters
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(-0.3).margin(tol));
    // REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(-0.3).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << -5.37070541971519e-06,  5.58776046545164e-06, -8.62056987828971e-07,
		-5.32137209241095e-06, -4.87624482275744e-06,  3.51720726336518e-06,
		 1.80217660206107e-06,  2.74275008576624e-06,  7.02453977468668e-06;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));

}


TEST_CASE("Wadenback 3DV 2026 - two-sided") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.093972209335993,  0.0431909090094019),
        Eigen::Vector2d(-0.0969912533237537,  0.0171160938161422),
        Eigen::Vector2d(-0.123299663772998,  0.500859415499538),
        Eigen::Vector2d(-0.300610529524079,  0.477995136756032),
        Eigen::Vector2d( 0.143634675965547, 0.402159621931429)
    };

    std::vector<Eigen::Vector2d> p2 = {
        Eigen::Vector2d(-0.0150485985127741,   0.504312576320796),
        Eigen::Vector2d(-0.0324017385815406,   0.525321116499528),
        Eigen::Vector2d(0.335015888737872, 0.182133325842435),
        Eigen::Vector2d(0.448124752741167, 0.313117448200805),
        Eigen::Vector2d(0.0822796392160207,  0.0746306773026531)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::Wadenback3DV2026::get_double_sided(p1, p2, true);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 1);

    // Test distortion parameters
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(-0.3).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(-0.3).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << -5.37070541971519e-06,  5.58776046545164e-06, -8.62056987828971e-07,
		-5.32137209241095e-06, -4.87624482275744e-06,  3.51720726336518e-06,
		 1.80217660206107e-06,  2.74275008576624e-06,  7.02453977468668e-06;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));

}
