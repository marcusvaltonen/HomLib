#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_nakano_icpr_2025.hpp"
#include "get_valtonenornhag_icpr_2026.hpp"
#include "get_wadenback_3dv_2026.hpp"
#include "flipped_solver.h"
#include "posedata.hpp"
#include "problem_instance.hpp"
#include "generate_problem_instance.hpp"

namespace {

HomLib::ProblemInstance make_instance(HomLib::DistortionCase dc, int n_points, unsigned int seed) {
    HomLib::ProblemConfig config;
    config.distortion = dc;
    config.point_noise = 0.0;
    config.number_points = n_points;
    config.seed = seed;
    return HomLib::generate_problem_instance(config);
}

bool has_gt_solution(const std::vector<HomLib::PoseData>& poses, const HomLib::ProblemInstance& inst,
                     double hom_tol = 1e-6, double dist_tol = 1e-6) {
    for (const auto& p : poses) {
        if (inst.hom_error(p.homography) < hom_tol &&
            inst.dist_error(p.distortion_parameter, p.distortion_parameter2) < dist_tol)
            return true;
    }
    return false;
}

template <typename NativeSolver>
void verify_flipped_on_valid_samples() {
    auto inst = make_instance(HomLib::ONE_SIDED_LEFT, 5, 42);

    NativeSolver native;
    HomLib::FlippedOneSidedSolver<NativeSolver> flipped;

    std::vector<HomLib::PoseData> native_poses, flipped_poses;
    native.solve(inst.x1, inst.x2, &native_poses);
    flipped.solve(inst.x2, inst.x1, &flipped_poses);

    REQUIRE(!native_poses.empty());
    REQUIRE(has_gt_solution(native_poses, inst));

    REQUIRE(native_poses.size() == flipped_poses.size());

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

TEST_CASE("Flipped Fitzgibbon: left and right are consistent on valid samples") {
    verify_flipped_on_valid_samples<HomLib::FitzgibbonCVPR2001::SolverSingleSided>();
}

TEST_CASE("Flipped Nakano: left and right are consistent on valid samples") {
    verify_flipped_on_valid_samples<HomLib::NakanoICPR2025::SolverSingleSided>();
}

TEST_CASE("Flipped Wadenback: left and right are consistent on valid samples") {
    verify_flipped_on_valid_samples<HomLib::Wadenback3DV2026::SolverSingleSided>();
}

TEST_CASE("Flipped Affine: left and right are consistent on valid samples") {
    auto inst = make_instance(HomLib::ONE_SIDED_RIGHT, 2, 43);

    using Native = HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided;
    Native native;
    HomLib::FlippedAffineOneSidedSolver<Native> flipped;

    std::vector<Eigen::Matrix2d> A_inv(inst.A.size());
    for (size_t i = 0; i < inst.A.size(); ++i)
        A_inv[i] = inst.A[i].inverse();

    std::vector<HomLib::PoseData> native_poses, flipped_poses;
    native.solve(inst.x1, inst.x2, inst.A, &native_poses);
    flipped.solve(inst.x2, inst.x1, A_inv, &flipped_poses);

    REQUIRE(!native_poses.empty());
    REQUIRE(has_gt_solution(native_poses, inst));

    REQUIRE(native_poses.size() == flipped_poses.size());

    for (size_t i = 0; i < native_poses.size(); ++i) {
        const auto& np = native_poses[i];
        const auto& fp = flipped_poses[i];

        REQUIRE(fp.distortion_parameter == Catch::Approx(np.distortion_parameter2).margin(1e-10));
        REQUIRE(fp.distortion_parameter2 == Catch::Approx(np.distortion_parameter).margin(1e-10));

        Eigen::Matrix3d H_inv = np.homography.inverse();
        REQUIRE(fp.homography.isApprox(H_inv, 1e-10));
    }
}

TEST_CASE("Flipped Orientation: left and right are consistent on valid samples") {
    auto inst = make_instance(HomLib::ONE_SIDED_RIGHT, 4, 44);

    using Native = HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided;
    Native native;
    HomLib::FlippedOrientationOneSidedSolver<Native> flipped;

    std::vector<Eigen::Vector2d> ori_swapped(inst.ori.size());
    for (size_t i = 0; i < inst.ori.size(); ++i)
        ori_swapped[i] = Eigen::Vector2d(inst.ori[i](1), inst.ori[i](0));

    std::vector<HomLib::PoseData> native_poses, flipped_poses;
    native.solve(inst.x1, inst.x2, inst.ori, &native_poses);
    flipped.solve(inst.x2, inst.x1, ori_swapped, &flipped_poses);

    REQUIRE(!native_poses.empty());
    REQUIRE(has_gt_solution(native_poses, inst));

    REQUIRE(native_poses.size() == flipped_poses.size());

    for (size_t i = 0; i < native_poses.size(); ++i) {
        const auto& np = native_poses[i];
        const auto& fp = flipped_poses[i];

        REQUIRE(fp.distortion_parameter == Catch::Approx(np.distortion_parameter2).margin(1e-10));
        REQUIRE(fp.distortion_parameter2 == Catch::Approx(np.distortion_parameter).margin(1e-10));

        Eigen::Matrix3d H_inv = np.homography.inverse();
        REQUIRE(fp.homography.isApprox(H_inv, 1e-10));
    }
}
