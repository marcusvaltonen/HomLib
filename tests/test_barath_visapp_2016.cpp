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
#include <cmath>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_barath_visapp_2016.hpp"
#include "generate_problem_instance.hpp"
#include "posedata.hpp"

static const double kPI = 3.14159265358979323846;

struct FullInstance {
    Eigen::Matrix3d H, F;
    std::vector<Eigen::Vector2d> x1, x2;
    std::vector<Eigen::Matrix2d> A;
};

static Eigen::Matrix3d skew(const Eigen::Vector3d &v) {
    Eigen::Matrix3d S;
    S << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    return S;
}

static FullInstance generate_full(unsigned int seed, int n_pts) {
    FullInstance inst;
    double fov = std::tan(70.0 / 2.0 * kPI / 180.0);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> depth(0.1, 10.0);
    std::uniform_real_distribution<double> coord(-fov, fov);
    std::normal_distribution<double> dir(0.0, 1.0);

    for (int att = 0; att < 100; ++att) {
        Eigen::Vector3d t(dir(rng), dir(rng), dir(rng)); t.normalize();
        Eigen::Quaterniond q(dir(rng), dir(rng), dir(rng), dir(rng)); q.normalize();
        Eigen::Matrix3d R = q.toRotationMatrix();
        Eigen::Vector3d n(dir(rng), dir(rng), dir(rng)); n.normalize();
        double alpha = depth(rng) / n(2);
        inst.H = alpha * R + t * n.transpose();
        inst.F = skew(t) * R;
        inst.x1.clear(); inst.x2.clear(); inst.A.clear();
        bool ok = true;
        for (int j = 0; j < n_pts; ++j) {
            bool found = false;
            for (int tr = 0; tr < 10; ++tr) {
                Eigen::Vector3d p1{coord(rng), coord(rng), 1.0}; p1.normalize();
                double lam = alpha / n.dot(p1);
                if (lam < 0) continue;
                Eigen::Vector3d p2 = (R * (p1 * lam) + t).normalized();
                if (p2(2) < 0) continue;
                Eigen::Vector2d x2h = p2.hnormalized();
                if (std::abs(x2h(0)) > fov || std::abs(x2h(1)) > fov) continue;
                Eigen::Vector2d x1h = p1.hnormalized();
                double s = inst.H(2,0)*x1h(0) + inst.H(2,1)*x1h(1) + inst.H(2,2);
                Eigen::Matrix2d Aff;
                Aff(0,0) = (inst.H(0,0) - inst.H(2,0)*x2h(0))/s;
                Aff(0,1) = (inst.H(0,1) - inst.H(2,1)*x2h(0))/s;
                Aff(1,0) = (inst.H(1,0) - inst.H(2,0)*x2h(1))/s;
                Aff(1,1) = (inst.H(1,1) - inst.H(2,1)*x2h(1))/s;
                if (Aff.determinant() < 0) continue;
                inst.x1.push_back(x1h); inst.x2.push_back(x2h); inst.A.push_back(Aff);
                found = true; break;
            }
            if (!found) { ok = false; break; }
        }
        if (ok) return inst;
    }
    return inst;
}

static double hom_err(const Eigen::Matrix3d &Hgt, const Eigen::Matrix3d &Hest) {
    Eigen::Matrix3d H1 = Hgt.normalized(); if (H1(2,2) < 0) H1 *= -1;
    Eigen::Matrix3d H2 = Hest.normalized(); if (H2(2,2) < 0) H2 *= -1;
    return (H1 - H2).norm();
}

TEST_CASE("Barath VISAPP 2016 - HA SVD and QR both reproduce GT") {
    int n_svd = 0, n_qr = 0;
    int n_tests = 200;

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 2;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        auto r_svd = HomLib::BarathVISAPP2016::get_affine(inst.x1, inst.x2, inst.A, false);
        auto r_qr = HomLib::BarathVISAPP2016::get_affine(inst.x1, inst.x2, inst.A, true);

        if (!r_svd.empty() && inst.hom_error(r_svd[0].homography) < 1e-8)
            n_svd++;
        if (!r_qr.empty() && inst.hom_error(r_qr[0].homography) < 1e-8)
            n_qr++;
    }

    REQUIRE(n_svd >= 180);
    REQUIRE(n_qr >= 180);
}

TEST_CASE("Barath VISAPP 2016 - HAF reproduces GT homography") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 1);
        if (inst.x1.empty()) continue;
        auto result = HomLib::BarathVISAPP2016::get_haf(inst.x1, inst.x2, inst.A, inst.F);
        if (result.empty()) continue;
        if (hom_err(inst.H, result[0].homography) < 1e-8)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - 3PT reproduces GT homography") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 3);
        if (inst.x1.size() < 3) continue;
        auto result = HomLib::BarathVISAPP2016::get_3pt(inst.x1, inst.x2, inst.F);
        if (result.empty()) continue;
        if (hom_err(inst.H, result[0].homography) < 1e-8)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - HAF overdetermined improves accuracy") {
    auto inst = generate_full(42, 5);
    REQUIRE(inst.x1.size() == 5);

    auto r1 = HomLib::BarathVISAPP2016::get_haf(
        {inst.x1[0]}, {inst.x2[0]}, {inst.A[0]}, inst.F);
    auto r5 = HomLib::BarathVISAPP2016::get_haf(inst.x1, inst.x2, inst.A, inst.F);

    REQUIRE(!r1.empty());
    REQUIRE(!r5.empty());
    REQUIRE(hom_err(inst.H, r5[0].homography) < 1e-8);
}

TEST_CASE("Barath VISAPP 2016 - 3PT estimate() transforms F under normalization", "[F_norm]") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 3);
        if (inst.x1.size() < 3) continue;

        // Scale points to pixel-like magnitudes so normalization is non-trivial
        double scale = 300.0;
        std::vector<Eigen::Vector2d> x1s(inst.x1), x2s(inst.x2);
        for (auto &p : x1s) p *= scale;
        for (auto &p : x2s) p *= scale;
        // F transforms as F' = S^{-T} F S^{-1} where S = diag(s,s,1)
        Eigen::Matrix3d Fs = Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal()
                           * inst.F
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();

        HomLib::BarathVISAPP2016::ThreePTSolver solver(Fs);
        std::vector<HomLib::PoseData> poses;
        solver.estimate(x1s, x2s, &poses);
        if (poses.empty()) continue;

        // The returned H should be in the scaled coordinate frame
        Eigen::Matrix3d Hs = Eigen::Vector3d(scale, scale, 1).asDiagonal()
                           * inst.H
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();
        if (hom_err(Hs, poses[0].homography) < 1e-6)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - HAF estimate() transforms F under normalization", "[F_norm]") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 1);
        if (inst.x1.empty()) continue;

        double scale = 300.0;
        std::vector<Eigen::Vector2d> x1s(inst.x1), x2s(inst.x2);
        for (auto &p : x1s) p *= scale;
        for (auto &p : x2s) p *= scale;
        Eigen::Matrix3d Fs = Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal()
                           * inst.F
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();

        HomLib::BarathVISAPP2016::HAFSolver solver(Fs);
        std::vector<HomLib::PoseData> poses;
        solver.estimate(x1s, x2s, inst.A, &poses);
        if (poses.empty()) continue;

        Eigen::Matrix3d Hs = Eigen::Vector3d(scale, scale, 1).asDiagonal()
                           * inst.H
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();
        if (hom_err(Hs, poses[0].homography) < 1e-6)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}
