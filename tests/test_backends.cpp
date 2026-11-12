#define _USE_MATH_DEFINES

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
#include <algorithm>
#include <cmath>
#include <limits>
#include <cfenv>
#include <catch2/catch_test_macros.hpp>

#include "posedata.hpp"
#include "radial.hpp"
#include "problem_instance.hpp"
#include "generate_problem_instance.hpp"
#include "get_wadenback_3dv_2026.hpp"
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_kukelova_cvpr_2015.hpp"
#include "get_valtonenornhag_icpr_2026.hpp"
#include "pose_estimator.h"

#include "lomsac_runner.h"
#include "solver_utils.h"

#ifdef HOMLIB_HAS_GCRANSAC
#include "gcransac_runner.h"
#endif

#ifdef HOMLIB_HAS_MAGSAC
#include "magsac_runner.h"
#endif

#ifdef HOMLIB_HAS_SUPERANSAC
#include "superansac_runner.h"
#endif


namespace {

struct TestData {
    HomLib::ProblemInstance instance;
    std::vector<Eigen::Vector2d> x;
    std::vector<Eigen::Vector2d> y;
    int n_total;
    int n_outliers;
};

TestData generate_test_data(int seed, int n_points, int n_outliers, HomLib::DistortionCase distortion,
                           double noise = 1e-4) {
    HomLib::ProblemConfig config;
    config.distortion = distortion;
    config.point_noise = noise;
    config.number_points = n_points;
    config.seed = static_cast<unsigned int>(seed);

    HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

    std::vector<int> indices(n_points);
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 rng(seed);
    std::shuffle(indices.begin(), indices.end(), rng);

    std::normal_distribution<double> outlier_dist(0.0, 5.0);
    for (int i = 0; i < n_outliers; ++i) {
        int idx = indices[i];
        inst.x1[idx] = Eigen::Vector2d(outlier_dist(rng), outlier_dist(rng));
        inst.x2[idx] = Eigen::Vector2d(outlier_dist(rng), outlier_dist(rng));
    }

    TestData data;
    data.instance = inst;
    data.x = inst.x1;
    data.y = inst.x2;
    data.n_total = n_points;
    data.n_outliers = n_outliers;
    return data;
}

} // anonymous namespace


// ---------------------------------------------------------------------------
// The residual that every backend uses
// ---------------------------------------------------------------------------

namespace {

// Points where computing the residual used to give a NaN or divide by zero
struct ResidualCase {
    const char* what;
    Eigen::Matrix3d H;
    double k1, k2;
    Eigen::Vector2d x;
    bool outlier;
};

std::vector<ResidualCase> residual_cases() {
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d H_inf = I;   // H z at infinity where z(0) = 1
    H_inf.row(2) << 1.0, 0.0, -1.0;
    return {
        {"4 k2 |H z|^2 = 1.6 > 1: H z can't be distorted", I, 0.0, 0.1, {2.0, 0.0}, true},
        {"H z at infinity", H_inf, 0.0, 0.0, {1.0, 0.5}, true},
        {"1 + k1 |x|^2 = 0: x can't be undistorted", I, -0.25, 0.0, {2.0, 0.0}, true},
        {"x = H z = 0", I, 0.0, 0.1, {0.0, 0.0}, false},
    };
}

} // anonymous namespace

TEST_CASE("Residual - a point without a residual is an outlier") {
    HomLib::Wadenback3DV2026::SolverTwoSided solver;
    const Eigen::Vector2d y(0.3, 0.0);
    for (const auto& c : residual_cases()) {
        INFO(c.what);
        HomLib::PoseData pose;
        pose.homography = c.H;
        pose.distortion_parameter = c.k1;
        pose.distortion_parameter2 = c.k2;
        const double r2 = HomLib::detail::compute_squared_residual(solver, pose, c.x, y);
        if (c.outlier)
            CHECK(r2 == std::numeric_limits<double>::max());
        else
            CHECK(r2 == (y - c.x).squaredNorm());
    }

    // Where the distortion can be inverted, the residual is unchanged
    HomLib::PoseData pose;   // H = I, no distortion in the first image
    pose.distortion_parameter2 = 0.1;
    const Eigen::Vector2d x1(1.0, 0.0);
    Eigen::Vector2d y1;
    REQUIRE(HomLib::try_radialdistort(x1, pose.distortion_parameter2, y1));
    REQUIRE(HomLib::detail::compute_squared_residual(solver, pose, x1, y1) < 1e-24);
    REQUIRE(HomLib::detail::compute_squared_residual(solver, pose, x1, x1)
            == (y1 - x1).squaredNorm());
}

TEST_CASE("Residual - no invalid operation or division by zero") {
    // With floating-point traps enabled (feenableexcept), either would raise SIGFPE. The inputs
    // are read through volatile after clearing the flags, and the result stored to a volatile
    // before testing them, so that the compiler neither folds nor moves the computation.
    HomLib::Wadenback3DV2026::SolverTwoSided solver;
    const Eigen::Vector2d y(0.3, 0.0);
    for (const auto& c : residual_cases()) {
        INFO(c.what);
        volatile double k1 = c.k1, k2 = c.k2, x0 = c.x(0), x1 = c.x(1);
        std::feclearexcept(FE_ALL_EXCEPT);
        HomLib::PoseData pose;
        pose.homography = c.H;
        pose.distortion_parameter = k1;
        pose.distortion_parameter2 = k2;
        volatile double r2 = HomLib::detail::compute_squared_residual(
            solver, pose, Eigen::Vector2d(x0, x1), y);
        const int flags = std::fetestexcept(FE_INVALID | FE_DIVBYZERO);
        (void)r2;
        CHECK((flags & FE_INVALID) == 0);
        CHECK((flags & FE_DIVBYZERO) == 0);
    }
}

TEST_CASE("LOMSAC - outliers where the true distortion can't be inverted") {
    // Pincushion distortion in the second image (k2 > 0): its distort has no real solution beyond
    // r = 1 / (2 sqrt(k2)) = 1.118. The inliers lie in the unit disk, and H maps them within that radius;
    // the outliers' x lie at r in [1.5, 2], so that H x is beyond 1.118 for the true model and
    // every model close to it. Their residual used to be NaN, which made the MSAC score of these
    // models NaN, so LOMSAC never selected them.
    const double k2 = 0.2;
    const double c = std::cos(0.2), s = std::sin(0.2);
    Eigen::Matrix3d H;
    H << c, -s, 0.05,
         s,  c, -0.03,
         0.01, 0.02, 1.0;
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> angle(0.0, 2.0 * M_PI), unit(0.0, 1.0);
    std::vector<Eigen::Vector2d> x, y;
    const int n_inliers = 160, n_outliers = 40;
    for (int i = 0; i < n_inliers + n_outliers; ++i) {
        const double a = angle(rng);
        const double r = i < n_inliers ? std::sqrt(unit(rng)) : 1.5 + 0.5 * unit(rng);
        const Eigen::Vector2d xi = r * Eigen::Vector2d(std::cos(a), std::sin(a));
        Eigen::Vector2d yi;
        if (i < n_inliers) {
            REQUIRE(HomLib::try_radialdistort((H * xi.homogeneous()).hnormalized(), k2, yi));
        } else {
            const double b = angle(rng), rb = std::sqrt(unit(rng));
            yi = rb * Eigen::Vector2d(std::cos(b), std::sin(b));
        }
        x.push_back(xi);
        y.push_back(yi);
    }
    HomLib::Wadenback3DV2026::SolverSingleSided solver;
    for (int i = 0; i < n_inliers + n_outliers; ++i) {
        const Eigen::Vector2d Hx = (H * x[i].homogeneous()).hnormalized();
        REQUIRE((4.0 * k2 * Hx.squaredNorm() < 1.0) == (i < n_inliers));
    }

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 7;
    auto [pose, stats] = HomLib::run_lomsac_point(solver, x, y, {}, options);

    REQUIRE(stats.best_num_inliers == n_inliers);
    REQUIRE(std::abs(pose.distortion_parameter2 - k2) < 1e-6);
    Eigen::Matrix3d H1 = H / H(2, 2), H2 = pose.homography / pose.homography(2, 2);
    REQUIRE((H1 - H2).norm() < 1e-6);
}


// ---------------------------------------------------------------------------
// LOMSAC — always available
// ---------------------------------------------------------------------------

TEST_CASE("LOMSAC - point features") {
    auto data = generate_test_data(42, 500, 100, HomLib::ONE_SIDED_LEFT);
    HomLib::Wadenback3DV2026::SolverSingleSided solver;

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 42;

    auto [pose, stats] = HomLib::run_lomsac_point(solver, data.x, data.y, {}, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.9 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("LOMSAC - affine features") {
    auto data = generate_test_data(43, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided solver;

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 43;

    auto [pose, stats] = HomLib::run_lomsac_affine(
        solver, data.x, data.y, data.instance.A, {}, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.9 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("LOMSAC - orientation features") {
    auto data = generate_test_data(44, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided solver;

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 44;

    auto [pose, stats] = HomLib::run_lomsac_orientation(
        solver, data.x, data.y, data.instance.ori, {}, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.9 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("LOMSAC - overdetermined non-minimal") {
    // Higher noise so non-minimal strategy actually affects the result
    auto data = generate_test_data(42, 500, 100, HomLib::ONE_SIDED_LEFT, 1e-3);
    HomLib::FitzgibbonCVPR2001::SolverSingleSided solver;

    HomLib::NonMinimalSettings nm;
    nm.strategy = HomLib::NonMinimalStrategy::Overdetermined;

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.003 * 0.003;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 42;

    auto [pose, stats] = HomLib::run_lomsac_point(solver, data.x, data.y, nm, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.9 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("LOMSAC - overdetermined affine") {
    auto data = generate_test_data(43, 200, 40, HomLib::ONE_SIDED_RIGHT, 1e-3);
    HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided solver;

    HomLib::NonMinimalSettings nm;
    nm.strategy = HomLib::NonMinimalStrategy::Overdetermined;

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 43;

    auto [pose, stats] = HomLib::run_lomsac_affine(
        solver, data.x, data.y, data.instance.A, nm, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.8 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("LOMSAC - overdetermined orientation") {
    auto data = generate_test_data(44, 200, 40, HomLib::ONE_SIDED_RIGHT, 1e-3);
    HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided solver;

    HomLib::NonMinimalSettings nm;
    nm.strategy = HomLib::NonMinimalStrategy::Overdetermined;

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 44;

    auto [pose, stats] = HomLib::run_lomsac_orientation(
        solver, data.x, data.y, data.instance.ori, nm, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.8 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("LOMSAC - overdetermined fallback for non-OD solver") {
    // Kukelova doesn't support overdetermined; the dispatch should
    // fall back to Sample without crashing or asserting.
    auto data = generate_test_data(42, 500, 100, HomLib::TWO_SIDED_EQUAL, 1e-3);
    HomLib::KukelovaCVPR2015::SolverTwoSidedEqual solver;

    HomLib::NonMinimalSettings nm;
    nm.strategy = HomLib::NonMinimalStrategy::Overdetermined;

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 42;

    auto [pose, stats] = HomLib::run_lomsac_point(solver, data.x, data.y, nm, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.8 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}


TEST_CASE("LOMSAC - overdetermined beats sample on average (no refine)") {
    HomLib::FitzgibbonCVPR2001::SolverSingleSided solver;

    int od_wins = 0;
    const int n_seeds = 5;
    for (int seed = 0; seed < n_seeds; ++seed) {
        auto data = generate_test_data(seed, 500, 100, HomLib::ONE_SIDED_LEFT, 1e-3);

        ransac_lib::LORansacOptions options;
        options.squared_inlier_threshold_ = 0.003 * 0.003;
        options.max_num_iterations_ = 5000;
        options.random_seed_ = seed;

        HomLib::NonMinimalSettings nm_sample;
        nm_sample.strategy = HomLib::NonMinimalStrategy::Sample;
        nm_sample.refine = false;
        auto [pose_s, stats_s] = HomLib::run_lomsac_point(
            solver, data.x, data.y, nm_sample, options);

        HomLib::NonMinimalSettings nm_od;
        nm_od.strategy = HomLib::NonMinimalStrategy::Overdetermined;
        nm_od.refine = false;
        auto [pose_od, stats_od] = HomLib::run_lomsac_point(
            solver, data.x, data.y, nm_od, options);

        double err_sample = data.instance.hom_error(pose_s.homography);
        double err_od = data.instance.hom_error(pose_od.homography);
        if (err_od < err_sample)
            ++od_wins;
    }
    REQUIRE(od_wins >= 4);
}


// ---------------------------------------------------------------------------
// GC-RANSAC — optional
// ---------------------------------------------------------------------------

#ifdef HOMLIB_HAS_GCRANSAC

TEST_CASE("GC-RANSAC - point features") {
    auto data = generate_test_data(42, 500, 100, HomLib::ONE_SIDED_LEFT);
    HomLib::Wadenback3DV2026::SolverSingleSided solver;

    HomLib::GCRANSACSettings settings;
    settings.threshold = 0.005;
    settings.max_iterations = 5000;

    // GC-RANSAC uses std::random_device internally (no seed API), so retry
    // up to 3 times to guard against unlucky sampling draws.
    size_t expected_inliers = data.n_total - data.n_outliers;
    bool success = false;
    for (int attempt = 0; attempt < 3 && !success; ++attempt) {
        auto [pose, stats] = HomLib::run_gcransac_point(solver, data.x, data.y, {}, {}, settings);
        success = (stats.inliers.size() > 0.9 * expected_inliers) &&
                  (data.instance.hom_error(pose.homography) < 0.1);
    }
    REQUIRE(success);
}

TEST_CASE("GC-RANSAC - affine features") {
    auto data = generate_test_data(43, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided solver;

    HomLib::GCRANSACSettings settings;
    settings.threshold = 0.005;
    settings.max_iterations = 5000;

    auto [pose, stats] = HomLib::run_gcransac_affine(
        solver, data.x, data.y, data.instance.A, {}, {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("GC-RANSAC - orientation features") {
    auto data = generate_test_data(44, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided solver;

    HomLib::GCRANSACSettings settings;
    settings.threshold = 0.005;
    settings.max_iterations = 5000;

    auto [pose, stats] = HomLib::run_gcransac_orientation(
        solver, data.x, data.y, data.instance.ori, {}, {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

#endif // HOMLIB_HAS_GCRANSAC


// ---------------------------------------------------------------------------
// MAGSAC++ — optional
// ---------------------------------------------------------------------------

#ifdef HOMLIB_HAS_MAGSAC

TEST_CASE("MAGSAC++ - point features") {
    auto data = generate_test_data(42, 500, 100, HomLib::ONE_SIDED_LEFT);
    HomLib::Wadenback3DV2026::SolverSingleSided solver;

    HomLib::MAGSACSettings settings;
    settings.sigma_max = 0.01;
    settings.max_iterations = 5000;

    auto [pose, stats] = HomLib::run_magsac_point(solver, data.x, data.y, {}, {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("MAGSAC++ - affine features") {
    auto data = generate_test_data(43, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided solver;

    HomLib::MAGSACSettings settings;
    settings.sigma_max = 0.01;
    settings.max_iterations = 5000;

    auto [pose, stats] = HomLib::run_magsac_affine(
        solver, data.x, data.y, data.instance.A, {}, {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("MAGSAC++ - orientation features") {
    auto data = generate_test_data(44, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided solver;

    HomLib::MAGSACSettings settings;
    settings.sigma_max = 0.01;
    settings.max_iterations = 5000;

    auto [pose, stats] = HomLib::run_magsac_orientation(
        solver, data.x, data.y, data.instance.ori, {}, {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

#endif // HOMLIB_HAS_MAGSAC


// ---------------------------------------------------------------------------
// SupeRANSAC — optional
// ---------------------------------------------------------------------------

#ifdef HOMLIB_HAS_SUPERANSAC

TEST_CASE("SupeRANSAC - point features") {
    auto data = generate_test_data(42, 500, 100, HomLib::ONE_SIDED_LEFT);
    HomLib::Wadenback3DV2026::SolverSingleSided solver;

    superansac::RANSACSettings settings;
    settings.inlierThreshold = 0.005;
    settings.maxIterations = 5000;
    settings.sampler = superansac::samplers::SamplerType::Uniform;

    auto [pose, stats] = HomLib::run_superansac_point(solver, data.x, data.y, {},
        {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("SupeRANSAC - affine features") {
    auto data = generate_test_data(43, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided solver;

    superansac::RANSACSettings settings;
    settings.inlierThreshold = 0.005;
    settings.maxIterations = 5000;
    settings.sampler = superansac::samplers::SamplerType::Uniform;

    auto [pose, stats] = HomLib::run_superansac_affine(
        solver, data.x, data.y, data.instance.A, {}, {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("SupeRANSAC - orientation features") {
    auto data = generate_test_data(44, 200, 40, HomLib::ONE_SIDED_RIGHT);
    HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided solver;

    superansac::RANSACSettings settings;
    settings.inlierThreshold = 0.005;
    settings.maxIterations = 5000;
    settings.sampler = superansac::samplers::SamplerType::Uniform;

    auto [pose, stats] = HomLib::run_superansac_orientation(
        solver, data.x, data.y, data.instance.ori, {}, {}, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

TEST_CASE("SupeRANSAC - the uniform sampler draws only indices of points") {
    // SupeRANSAC's UniformRandomSampler, initialized on the data, can draw the index of the row
    // after the last point; create_sampler gives HomLib's
    const DataMatrix data = DataMatrix::Zero(5, 4);
    auto sampler = HomLib::detail::create_sampler<4>(superansac::samplers::SamplerType::Uniform);
    sampler->initialize(data);
    size_t sample[4];
    bool in_range = true;
    for (int i = 0; i < 1000 && in_range; i++) {
        REQUIRE(sampler->sample(data, 4, sample));
        for (size_t j : sample)
            in_range = in_range && j < 5;
    }
    REQUIRE(in_range);
}

#endif // HOMLIB_HAS_SUPERANSAC
