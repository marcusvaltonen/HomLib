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
#include <random>
#include <cmath>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_nakano_icpr_2025.hpp"
#include "get_wadenback_3dv_2026.hpp"
#include "radial.hpp"
#include "posedata.hpp"
#include "solver_utils.h"

static void generate_one_sided_left(
    const Eigen::Matrix3d& H, double k2, int n,
    std::vector<Eigen::Vector2d>* x, std::vector<Eigen::Vector2d>* y,
    unsigned seed = 42)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> coord(-0.8, 0.8);
    x->resize(n);
    y->resize(n);
    for (int i = 0; i < n; i++) {
        (*x)[i] = Eigen::Vector2d(coord(rng), coord(rng));
        Eigen::Vector3d tmp = H * (*x)[i].homogeneous();
        Eigen::Vector2d yu = tmp.hnormalized();
        Eigen::Vector2d yd;
        REQUIRE(HomLib::try_radialdistort(yu, k2, yd));
        (*y)[i] = yd;
    }
}

static void generate_two_sided_equal(
    const Eigen::Matrix3d& H, double k, int n,
    std::vector<Eigen::Vector2d>* x, std::vector<Eigen::Vector2d>* y,
    unsigned seed = 42)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> coord(-0.8, 0.8);
    x->resize(n);
    y->resize(n);
    for (int i = 0; i < n; i++) {
        Eigen::Vector2d xu(coord(rng), coord(rng));
        Eigen::Vector2d xd, yd;
        REQUIRE(HomLib::try_radialdistort(xu, k, xd));
        (*x)[i] = xd;
        Eigen::Vector3d tmp = H * xu.homogeneous();
        Eigen::Vector2d yu = tmp.hnormalized();
        REQUIRE(HomLib::try_radialdistort(yu, k, yd));
        (*y)[i] = yd;
    }
}

// ============================================================================
// Fitzgibbon one-sided: get_single with N > 5
// ============================================================================

TEST_CASE("Fitzgibbon single - overdetermined N > 5") {
    Eigen::Matrix3d H;
    H << 0.9, 0.1, 0.05,
        -0.1, 0.85, 0.03,
         0.02, -0.01, 1.0;
    double k2 = -0.15;

    std::vector<Eigen::Vector2d> x, y;
    generate_one_sided_left(H, k2, 20, &x, &y);

    auto output = HomLib::FitzgibbonCVPR2001::get_single(x, y);
    REQUIRE(output.size() == 1);
    REQUIRE(output[0].distortion_parameter2 == Catch::Approx(k2).margin(1e-3));

    Eigen::Matrix3d H_est = output[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_gt = H / H(2, 2);
    double max_reproj = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        Eigen::Vector2d yu_est = (H_est * x[i].homogeneous()).hnormalized();
        Eigen::Vector2d yu_gt = (H_gt * x[i].homogeneous()).hnormalized();
        max_reproj = std::max(max_reproj, (yu_est - yu_gt).norm());
    }
    REQUIRE(max_reproj < 1e-3);
}

TEST_CASE("Fitzgibbon single - fewer than 5 returns empty") {
    std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d(0.1, 0.2), Eigen::Vector2d(0.3, 0.4),
        Eigen::Vector2d(0.5, 0.6), Eigen::Vector2d(0.7, 0.8)
    };
    REQUIRE(HomLib::FitzgibbonCVPR2001::get_single(x, x).empty());
}

// ============================================================================
// Fitzgibbon two-sided equal: get with N > 5
// ============================================================================

TEST_CASE("Fitzgibbon two-sided equal - overdetermined N > 5") {
    Eigen::Matrix3d H;
    H << 0.9, 0.1, 0.05,
        -0.1, 0.85, 0.03,
         0.02, -0.01, 1.0;
    double k = -0.15;

    std::vector<Eigen::Vector2d> x, y;
    generate_two_sided_equal(H, k, 20, &x, &y);

    auto output = HomLib::FitzgibbonCVPR2001::get(x, y);
    REQUIRE(output.size() == 1);
    REQUIRE(output[0].distortion_parameter == Catch::Approx(k).margin(1e-3));
    REQUIRE(output[0].distortion_parameter2 == Catch::Approx(k).margin(1e-3));
}

// ============================================================================
// Nakano: get with N > 5
// ============================================================================

TEST_CASE("Nakano - overdetermined N > 5") {
    Eigen::Matrix3d H;
    H << 0.9, 0.1, 0.05,
        -0.1, 0.85, 0.03,
         0.02, -0.01, 1.0;
    double k2 = -0.15;

    std::vector<Eigen::Vector2d> x, y;
    generate_one_sided_left(H, k2, 20, &x, &y);

    auto output = HomLib::NakanoICPR2025::get(x, y, true);
    REQUIRE(output.size() == 1);
    REQUIRE(output[0].distortion_parameter2 == Catch::Approx(k2).margin(1e-3));

    Eigen::Matrix3d H_est = output[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_gt = H / H(2, 2);
    double max_reproj = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        Eigen::Vector2d yu_est = (H_est * x[i].homogeneous()).hnormalized();
        Eigen::Vector2d yu_gt = (H_gt * x[i].homogeneous()).hnormalized();
        max_reproj = std::max(max_reproj, (yu_est - yu_gt).norm());
    }
    REQUIRE(max_reproj < 1e-3);
}

TEST_CASE("Nakano - fewer than 5 returns empty") {
    std::vector<Eigen::Vector2d> x = {
        Eigen::Vector2d(0.1, 0.2), Eigen::Vector2d(0.3, 0.4),
        Eigen::Vector2d(0.5, 0.6), Eigen::Vector2d(0.7, 0.8)
    };
    REQUIRE(HomLib::NakanoICPR2025::get(x, x, false).empty());
}

// ============================================================================
// Non-minimal dispatch: Wadenback minimal + Nakano non-minimal
// ============================================================================

TEST_CASE("non_minimal_dispatch_point with Nakano as non-minimal solver") {
    Eigen::Matrix3d H;
    H << 0.9, 0.1, 0.05,
        -0.1, 0.85, 0.03,
         0.02, -0.01, 1.0;
    double k2 = -0.15;

    std::vector<Eigen::Vector2d> x, y;
    generate_one_sided_left(H, k2, 20, &x, &y);

    HomLib::Wadenback3DV2026::SolverSingleSided minimal_solver;

    HomLib::NakanoICPR2025::SolverSingleSided nakano;
    HomLib::NonMinimalEstimateFn nm_fn = [&nakano](
        const std::vector<Eigen::Vector2d>& xv,
        const std::vector<Eigen::Vector2d>& yv,
        std::vector<HomLib::PoseData>* poses) -> int {
        nakano.estimate(xv, yv, poses);
        return static_cast<int>(poses->size());
    };

    auto est_fn = [](const HomLib::Wadenback3DV2026::SolverSingleSided& s,
                     const std::vector<Eigen::Vector2d>& xv,
                     const std::vector<Eigen::Vector2d>& yv,
                     const std::vector<size_t>& idx, size_t ss,
                     std::vector<HomLib::PoseData>* poses) {
        std::vector<Eigen::Vector2d> xs(ss), ys(ss);
        for (size_t i = 0; i < ss; ++i) { xs[i] = xv[idx[i]]; ys[i] = yv[idx[i]]; }
        s.estimate(xs, ys, poses);
    };

    HomLib::PoseData best;
    bool found = HomLib::detail::non_minimal_dispatch_point(
        HomLib::NonMinimalStrategy::Overdetermined,
        minimal_solver, x, y, est_fn, best, -1, nm_fn);

    REQUIRE(found);
    REQUIRE(best.distortion_parameter2 == Catch::Approx(k2).margin(1e-2));
}

TEST_CASE("non_minimal_dispatch_point falls back to Sample without fn") {
    Eigen::Matrix3d H;
    H << 0.9, 0.1, 0.05,
        -0.1, 0.85, 0.03,
         0.02, -0.01, 1.0;
    double k2 = -0.15;

    std::vector<Eigen::Vector2d> x, y;
    generate_one_sided_left(H, k2, 20, &x, &y);

    HomLib::Wadenback3DV2026::SolverSingleSided solver;

    auto est_fn = [](const HomLib::Wadenback3DV2026::SolverSingleSided& s,
                     const std::vector<Eigen::Vector2d>& xv,
                     const std::vector<Eigen::Vector2d>& yv,
                     const std::vector<size_t>& idx, size_t ss,
                     std::vector<HomLib::PoseData>* poses) {
        std::vector<Eigen::Vector2d> xs(ss), ys(ss);
        for (size_t i = 0; i < ss; ++i) { xs[i] = xv[idx[i]]; ys[i] = yv[idx[i]]; }
        s.estimate(xs, ys, poses);
    };

    HomLib::PoseData best;
    bool found = HomLib::detail::non_minimal_dispatch_point(
        HomLib::NonMinimalStrategy::Overdetermined,
        solver, x, y, est_fn, best);

    REQUIRE(found);
}
