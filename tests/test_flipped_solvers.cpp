#include <Eigen/Dense>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_nakano_icpr_2025.hpp"
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
