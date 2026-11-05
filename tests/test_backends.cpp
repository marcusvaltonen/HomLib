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
#include <catch2/catch_test_macros.hpp>

#include "posedata.hpp"
#include "problem_instance.hpp"
#include "generate_problem_instance.hpp"
#include "get_wadenback_3dv_2026.hpp"
#include "get_valtonenornhag_icpr_2026.hpp"

#include "lomsac_runner.h"

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

TestData generate_test_data(int seed, int n_points, int n_outliers, HomLib::DistortionCase distortion) {
    HomLib::ProblemConfig config;
    config.distortion = distortion;
    config.point_noise = 0.0001;
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
// LOMSAC — always available
// ---------------------------------------------------------------------------

TEST_CASE("LOMSAC - point features") {
    auto data = generate_test_data(42, 500, 100, HomLib::ONE_SIDED_LEFT);
    HomLib::Wadenback3DV2026::SolverSingleSided solver(true);

    ransac_lib::LORansacOptions options;
    options.squared_inlier_threshold_ = 0.005 * 0.005;
    options.max_num_iterations_ = 5000;
    options.random_seed_ = 42;

    auto [pose, stats] = HomLib::run_lomsac_point(solver, data.x, data.y, options);

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
        solver, data.x, data.y, data.instance.A, options);

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
        solver, data.x, data.y, data.instance.ori, options);

    double expected_inlier_ratio = 1.0 - static_cast<double>(data.n_outliers) / data.n_total;
    REQUIRE(stats.inlier_ratio > 0.9 * expected_inlier_ratio);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}


// ---------------------------------------------------------------------------
// GC-RANSAC — optional
// ---------------------------------------------------------------------------

#ifdef HOMLIB_HAS_GCRANSAC

TEST_CASE("GC-RANSAC - point features") {
    auto data = generate_test_data(42, 500, 100, HomLib::ONE_SIDED_LEFT);
    HomLib::Wadenback3DV2026::SolverSingleSided solver(true);

    HomLib::GCRANSACSettings settings;
    settings.threshold = 0.005;
    settings.max_iterations = 5000;

    // GC-RANSAC uses std::random_device internally (no seed API), so retry
    // up to 3 times to guard against unlucky sampling draws.
    size_t expected_inliers = data.n_total - data.n_outliers;
    bool success = false;
    for (int attempt = 0; attempt < 3 && !success; ++attempt) {
        auto [pose, stats] = HomLib::run_gcransac_point(solver, data.x, data.y, settings);
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
        solver, data.x, data.y, data.instance.A, settings);

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
        solver, data.x, data.y, data.instance.ori, settings);

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
    HomLib::Wadenback3DV2026::SolverSingleSided solver(true);

    HomLib::MAGSACSettings settings;
    settings.sigma_max = 0.01;
    settings.max_iterations = 5000;

    auto [pose, stats] = HomLib::run_magsac_point(solver, data.x, data.y, settings);

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
        solver, data.x, data.y, data.instance.A, settings);

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
        solver, data.x, data.y, data.instance.ori, settings);

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
    HomLib::Wadenback3DV2026::SolverSingleSided solver(true);

    superansac::RANSACSettings settings;
    settings.inlierThreshold = 0.005;
    settings.maxIterations = 5000;
    settings.sampler = superansac::samplers::SamplerType::Uniform;

    auto [pose, stats] = HomLib::run_superansac_point(solver, data.x, data.y, settings);

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
        solver, data.x, data.y, data.instance.A, settings);

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
        solver, data.x, data.y, data.instance.ori, settings);

    size_t expected_inliers = data.n_total - data.n_outliers;
    REQUIRE(stats.inliers.size() > 0.9 * expected_inliers);
    REQUIRE(data.instance.hom_error(pose.homography) < 0.1);
}

#endif // HOMLIB_HAS_SUPERANSAC
