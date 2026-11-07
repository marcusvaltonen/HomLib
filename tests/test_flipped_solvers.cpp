#include <Eigen/Dense>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_nakano_icpr_2025.hpp"
#include "get_valtonenornhag_icpr_2026.hpp"
#include "get_wadenback_3dv_2026.hpp"
#include "flipped_solver.h"
#include "posedata.hpp"

namespace {

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

template <typename NativeSolver>
void verify_flipped_consistency() {
    NativeSolver native;
    HomLib::FlippedOneSidedSolver<NativeSolver> flipped;

    std::vector<HomLib::PoseData> native_poses, flipped_poses;
    native.solve(p1, p2, &native_poses);
    flipped.solve(p2, p1, &flipped_poses);

    REQUIRE(native_poses.size() == flipped_poses.size());
    REQUIRE(!native_poses.empty());

    for (size_t i = 0; i < native_poses.size(); ++i) {
        const auto& np = native_poses[i];
        const auto& fp = flipped_poses[i];

        REQUIRE(fp.distortion_parameter == Catch::Approx(np.distortion_parameter2).margin(1e-10));
        REQUIRE(fp.distortion_parameter2 == Catch::Approx(np.distortion_parameter).margin(1e-10));

        Eigen::Matrix3d H_inv = np.homography.inverse();
        REQUIRE(fp.homography.isApprox(H_inv, 1e-10));
    }
}

} // anonymous namespace

TEST_CASE("Flipped Fitzgibbon: left and right are consistent") {
    verify_flipped_consistency<HomLib::FitzgibbonCVPR2001::SolverSingleSided>();
}

TEST_CASE("Flipped Nakano: left and right are consistent") {
    verify_flipped_consistency<HomLib::NakanoICPR2025::SolverSingleSided>();
}

TEST_CASE("Flipped Wadenback: left and right are consistent") {
    verify_flipped_consistency<HomLib::Wadenback3DV2026::SolverSingleSided>();
}

TEST_CASE("Flipped Affine: left and right are consistent") {
    using Native = HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided;
    Native native;
    HomLib::FlippedAffineOneSidedSolver<Native> flipped;

    std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d(-0.304372059145248, -0.174125682172663),
        Eigen::Vector2d(-0.308082306270758,  0.100342577802908)
    };
    std::vector<Eigen::Vector2d> y = {
        Eigen::Vector2d(0.596062339231747, 0.449851750693698),
        Eigen::Vector2d(0.284590553127076, 0.20125992440591)
    };
    std::vector<Eigen::Matrix2d> A = {
        (Eigen::Matrix2d() << -0.742737063351745, -1.25720051435581,
                               1.29931162730836, -0.953395533218031).finished(),
        (Eigen::Matrix2d() << -0.763314941736403, -1.05311577674003,
                               1.13911319830585, -0.839438880618645).finished()
    };

    // From the left-side perspective, affine frames are inverted
    std::vector<Eigen::Matrix2d> A_inv(A.size());
    for (size_t i = 0; i < A.size(); ++i)
        A_inv[i] = A[i].inverse();

    std::vector<HomLib::PoseData> native_poses, flipped_poses;
    native.solve(x, y, A, &native_poses);
    flipped.solve(y, x, A_inv, &flipped_poses);

    REQUIRE(native_poses.size() == flipped_poses.size());
    REQUIRE(!native_poses.empty());

    for (size_t i = 0; i < native_poses.size(); ++i) {
        const auto& np = native_poses[i];
        const auto& fp = flipped_poses[i];

        REQUIRE(fp.distortion_parameter == Catch::Approx(np.distortion_parameter2).margin(1e-10));
        REQUIRE(fp.distortion_parameter2 == Catch::Approx(np.distortion_parameter).margin(1e-10));

        Eigen::Matrix3d H_inv = np.homography.inverse();
        REQUIRE(fp.homography.isApprox(H_inv, 1e-10));
    }
}

TEST_CASE("Flipped Orientation: left and right are consistent") {
    using Native = HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided;
    Native native;
    HomLib::FlippedOrientationOneSidedSolver<Native> flipped;

    std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d(0.611733207361777, -0.100815813483642),
        Eigen::Vector2d(0.506500688761625,  0.22654451870654),
        Eigen::Vector2d(0.618492431042353,  0.490634389234033),
        Eigen::Vector2d(0.446207436531052, -0.221543454776215)
    };
    std::vector<Eigen::Vector2d> y = {
        Eigen::Vector2d( 0.152378076324475, 0.471173894309057),
        Eigen::Vector2d( 0.35249479911981,  0.604057206211103),
        Eigen::Vector2d( 0.576429540272504, 0.546190419126011),
        Eigen::Vector2d(-0.0960825885072252, 0.57183541266364)
    };
    std::vector<Eigen::Vector2d> ori = {
        Eigen::Vector2d(2.82616873462069, 2.11084467936878),
        Eigen::Vector2d(2.76670870212448, 1.5791745485192),
        Eigen::Vector2d(3.03912230291483, 1.40242855410085),
        Eigen::Vector2d(2.54398570702913, 2.11604504644583)
    };

    // From the left-side perspective, orientations are swapped
    std::vector<Eigen::Vector2d> ori_swapped(ori.size());
    for (size_t i = 0; i < ori.size(); ++i)
        ori_swapped[i] = Eigen::Vector2d(ori[i](1), ori[i](0));

    std::vector<HomLib::PoseData> native_poses, flipped_poses;
    native.solve(x, y, ori, &native_poses);
    flipped.solve(y, x, ori_swapped, &flipped_poses);

    REQUIRE(native_poses.size() == flipped_poses.size());
    REQUIRE(!native_poses.empty());

    for (size_t i = 0; i < native_poses.size(); ++i) {
        const auto& np = native_poses[i];
        const auto& fp = flipped_poses[i];

        REQUIRE(fp.distortion_parameter == Catch::Approx(np.distortion_parameter2).margin(1e-10));
        REQUIRE(fp.distortion_parameter2 == Catch::Approx(np.distortion_parameter).margin(1e-10));

        Eigen::Matrix3d H_inv = np.homography.inverse();
        REQUIRE(fp.homography.isApprox(H_inv, 1e-10));
    }
}
