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
#include "wadenback_common.hpp"


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
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(0.0).margin(tol));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(5.11091).margin(tol));
    REQUIRE(posedata[1].distortion_parameter == Catch::Approx(0.0).margin(tol));
    REQUIRE(posedata[1].distortion_parameter2 == Catch::Approx(2.73767).margin(tol));
    REQUIRE(posedata[2].distortion_parameter == Catch::Approx(0.0).margin(tol));
    REQUIRE(posedata[2].distortion_parameter2 == Catch::Approx(-0.3).margin(tol));

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
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(-0.3).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << -5.37070541971519e-06,  5.58776046545164e-06, -8.62056987828971e-07,
		-5.32137209241095e-06, -4.87624482275744e-06,  3.51720726336518e-06,
		 1.80217660206107e-06,  2.74275008576624e-06,  7.02453977468668e-06;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));

}


namespace {
// The points of the two-sided test (k1 = k2 = -0.3)
const std::vector<Eigen::Vector2d> kTwoSidedX = {
    Eigen::Vector2d(-0.093972209335993,  0.0431909090094019),
    Eigen::Vector2d(-0.0969912533237537,  0.0171160938161422),
    Eigen::Vector2d(-0.123299663772998,  0.500859415499538),
    Eigen::Vector2d(-0.300610529524079,  0.477995136756032),
    Eigen::Vector2d( 0.143634675965547, 0.402159621931429)
};
const std::vector<Eigen::Vector2d> kTwoSidedY = {
    Eigen::Vector2d(-0.0150485985127741,   0.504312576320796),
    Eigen::Vector2d(-0.0324017385815406,   0.525321116499528),
    Eigen::Vector2d(0.335015888737872, 0.182133325842435),
    Eigen::Vector2d(0.448124752741167, 0.313117448200805),
    Eigen::Vector2d(0.0822796392160207,  0.0746306773026531)
};
bool HasSolution(const std::vector<HomLib::PoseData> &poses, double k1, double k2, double tol) {
    for (const HomLib::PoseData &p : poses)
        if (std::abs(p.distortion_parameter - k1) < tol && std::abs(p.distortion_parameter2 - k2) < tol)
            return true;
    return false;
}
}  // namespace

TEST_CASE("Wadenback 3DV 2026 - two-sided returns all the real roots") {
    // The solver doesn't check the range of the distortion parameters (the estimator does, see the
    // next test): its roots include implausible ones
    std::vector<HomLib::PoseData> all = HomLib::Wadenback3DV2026::get_double_sided(kTwoSidedX, kTwoSidedY, false);
    REQUIRE(HasSolution(all, -0.3, -0.3, 1e-9));
    bool outside = false;
    for (const HomLib::PoseData &p : all)
        outside = outside || p.distortion_parameter < -1.0 || p.distortion_parameter > 0.1;
    REQUIRE(outside);
}

TEST_CASE("Wadenback 3DV 2026 - two-sided estimator with a range in the input units") {
    // The test problem in other units: the points scaled by s, so k scaled by 1/s^2
    const double s = 1000.0, k = -0.3 / (s * s);
    std::vector<Eigen::Vector2d> x, y;
    for (int i = 0; i < 5; i++) {
        x.push_back(s * kTwoSidedX[i]);
        y.push_back(s * kTwoSidedY[i]);
    }
    HomLib::Wadenback3DV2026::SolverTwoSided solver;
    std::vector<HomLib::PoseData> poses, unfiltered;

    // Default: on, in the normalized coordinates, whatever the units
    REQUIRE(solver.use_k_range);
    REQUIRE_FALSE(solver.k_range_in_input_units);
    solver.estimate(x, y, &poses);
    REQUIRE(HasSolution(poses, k, k, 1e-7 / (s * s)));
    solver.use_k_range = false;
    solver.estimate(x, y, &unfiltered);
    REQUIRE(unfiltered.size() > poses.size());
    solver.use_k_range = true;

    // Range in the input units, including or excluding k
    solver.k_range_in_input_units = true;
    solver.k_min = -0.5 / (s * s);
    solver.k_max = 0.0;
    solver.estimate(x, y, &poses);
    REQUIRE(HasSolution(poses, k, k, 1e-7 / (s * s)));
    for (const HomLib::PoseData &p : poses) {
        REQUIRE(p.distortion_parameter >= solver.k_min);
        REQUIRE(p.distortion_parameter <= solver.k_max);
        REQUIRE(p.distortion_parameter2 >= solver.k_min);
        REQUIRE(p.distortion_parameter2 <= solver.k_max);
    }
    solver.k_min = -0.2 / (s * s);
    solver.estimate(x, y, &poses);
    REQUIRE_FALSE(HasSolution(poses, k, k, 1e-7 / (s * s)));
}

TEST_CASE("Wadenback 3DV 2026 - one-sided estimator with a range check") {
    // The points of the one-sided test (k2 = -0.3; k1 = 0, not estimated)
    const std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d( 0.242104535827652, -0.296297354299645),
        Eigen::Vector2d( 0.319776775526576, -0.114580331777873),
        Eigen::Vector2d(-0.122936754370192, -0.532655969828257),
        Eigen::Vector2d( 0.303141152380151, -0.295669384794041),
        Eigen::Vector2d( 0.446905125282507, -0.528127859697035)
    };
    const std::vector<Eigen::Vector2d> y = {
        Eigen::Vector2d(-0.562787932824692, -0.0769625591786919),
        Eigen::Vector2d(-0.612813481767545, -0.0311462459576149),
        Eigen::Vector2d(-0.577247052139350, -0.252034297773712),
        Eigen::Vector2d(-0.544031097770073, -0.0440402280844913),
        Eigen::Vector2d(-0.424602042497788,  0.0211895908990797)
    };
    HomLib::Wadenback3DV2026::SolverSingleSided solver;
    solver.use_k_range = true;
    solver.k_range_in_input_units = true;
    solver.k_min = -0.5;   // excludes 0: the parameter that the solver doesn't estimate passes anyway
    solver.k_max = -0.1;
    std::vector<HomLib::PoseData> poses;
    solver.estimate(x, y, &poses);
    bool found = false;
    for (const HomLib::PoseData &p : poses)
        found = found || (p.distortion_parameter == 0 && std::abs(p.distortion_parameter2 + 0.3) < 1e-7);
    REQUIRE(found);
    solver.k_max = -0.4;   // excludes k2
    solver.estimate(x, y, &poses);
    REQUIRE(poses.empty());
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

    // Test that the true solution is among the solutions
    size_t i = 0;
    while (i < posedata.size() && std::abs(posedata[i].distortion_parameter + 0.3) > 1e-6)
        i++;
    REQUIRE(i < posedata.size());

    // Test distortion parameters
    REQUIRE(posedata[i].distortion_parameter == Catch::Approx(-0.3).margin(tol));
    REQUIRE(posedata[i].distortion_parameter2 == Catch::Approx(-0.3).margin(tol));

    // Test homographies
    tol = 1e-7;
    Eigen::Matrix3d expected;

    expected << -5.37070541971519e-06,  5.58776046545164e-06, -8.62056987828971e-07,
		-5.32137209241095e-06, -4.87624482275744e-06,  3.51720726336518e-06,
		 1.80217660206107e-06,  2.74275008576624e-06,  7.02453977468668e-06;
    REQUIRE(posedata[i].homography.isApprox(expected, tol));

}

namespace {
// [x{1} x{2} x{3}] with the division model w_i = 1 + k*r_i^2
Eigen::Matrix3d DistortedTriple(const std::vector<Eigen::Vector2d> &x, double k) {
    Eigen::Matrix3d X;
    for (int i = 0; i < 3; i++)
        X.col(i) << x[i](0), x[i](1), 1 + k * x[i].squaredNorm();
    return X;
}
// g = adj([x{1} x{2} x{3}]) * x{4}: the rows of the adjugate are x{2} x x{3}, x{3} x x{1}, x{1} x x{2}
Eigen::Vector3d AdjugateTimesFourth(const std::vector<Eigen::Vector2d> &x, double k) {
    const Eigen::Matrix3d X = DistortedTriple(x, k);
    Eigen::Matrix3d adj;
    adj.row(0) = X.col(1).cross(X.col(2)).transpose();
    adj.row(1) = X.col(2).cross(X.col(0)).transpose();
    adj.row(2) = X.col(0).cross(X.col(1)).transpose();
    const Eigen::Vector3d x4(x[3](0), x[3](1), 1 + k * x[3].squaredNorm());
    return adj * x4;
}
}  // namespace

TEST_CASE("Wadenback 3DV 2026 - residuals of the false-root checks") {
    const std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d(-0.093972209335993,  0.0431909090094019),
        Eigen::Vector2d(-0.0969912533237537,  0.0171160938161422),
        Eigen::Vector2d(-0.123299663772998,  0.500859415499538),
        Eigen::Vector2d(-0.300610529524079,  0.477995136756032),
        Eigen::Vector2d( 0.143634675965547, 0.402159621931429)
    };

    for (double k : {0.0, -0.05, -0.3, -0.7, 0.1}) {
        const Eigen::Matrix3d X = DistortedTriple(x, k);
        // compute_residual1 = |det([x{1} x{2} x{3}])| / product of the row norms
        const double r1 = std::abs(X.determinant()) / (X.row(0).norm() * X.row(1).norm() * X.row(2).norm());
        REQUIRE(HomLib::Wadenback3DV2026::compute_residual1(x, k) == Catch::Approx(r1).epsilon(1e-12));

        // compute_residual2 = min_i |g_i| / |g|
        const Eigen::Vector3d g = AdjugateTimesFourth(x, k);
        const double r2 = g.cwiseAbs().minCoeff() / g.norm();
        REQUIRE(HomLib::Wadenback3DV2026::compute_residual2(x, k) == Catch::Approx(r2).epsilon(1e-12));

        // Both are scale-free: the same for the points scaled by s and the distortion by 1/s^2
        const double s = 37.0;
        std::vector<Eigen::Vector2d> xs;
        for (const Eigen::Vector2d &p : x)
            xs.push_back(s * p);
        REQUIRE(HomLib::Wadenback3DV2026::compute_residual1(xs, k / (s * s)) == Catch::Approx(r1).epsilon(1e-10));
        REQUIRE(HomLib::Wadenback3DV2026::compute_residual2(xs, k / (s * s)) == Catch::Approx(r2).epsilon(1e-10));
    }

    // compute_residual2 vanishes where a component of g does (as at the false roots it detects,
    // where the same component vanishes in both images): g1 is affine in k
    const double g1_0 = AdjugateTimesFourth(x, 0.0)(0);
    const double g1_1 = AdjugateTimesFourth(x, 1.0)(0) - g1_0;
    const double k_false = -g1_0 / g1_1;
    REQUIRE(HomLib::Wadenback3DV2026::compute_residual2(x, k_false) < 1e-12);
}

TEST_CASE("Wadenback 3DV 2026 - two-sided keeps true solutions with small residuals") {
    // A noise-free problem (k1 = -0.165, k2 = -0.161) whose true solution was discarded by the
    // former absolute threshold on the residuals
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.11055214114931514, -0.5097570468713325),
        Eigen::Vector2d( 0.15869951293104498, -0.16205554084328394),
        Eigen::Vector2d(-0.13626033417641306, -0.48660931896153381),
        Eigen::Vector2d(-0.2301040297556369, -0.58014902645235755),
        Eigen::Vector2d(-0.3430612813896049, -0.6381657633318355)
    };
    std::vector<Eigen::Vector2d> p2 = {
        Eigen::Vector2d(-0.06561068159991755, 0.1735283567957854),
        Eigen::Vector2d( 0.30197423189591727, 0.46978937888994693),
        Eigen::Vector2d(-0.088204632690722554, 0.19930711144479341),
        Eigen::Vector2d(-0.19396490905569072, 0.12714507402566574),
        Eigen::Vector2d(-0.30856906316803173, 0.089876754041425103)
    };
    Eigen::Matrix3d expected;
    expected <<  8.738436073143756,   1.3439137301342539, 1.1335492294046015,
                -1.6550398019981574,  6.6971064929211241, 4.9572217965732133,
                -0.29241123402938313, -4.501701870435423, 6.5856346279614115;

    std::vector<HomLib::PoseData> posedata = HomLib::Wadenback3DV2026::get_double_sided(p1, p2, true);

    bool found = false;
    for (const HomLib::PoseData &pd : posedata) {
        Eigen::Matrix3d H = pd.homography / pd.homography(2, 2) * expected(2, 2);
        if (std::abs(pd.distortion_parameter + 0.16512790385947976) < 1e-9 &&
            std::abs(pd.distortion_parameter2 + 0.16057618398228576) < 1e-9 && H.isApprox(expected, 1e-7))
            found = true;
    }
    REQUIRE(found);
}

