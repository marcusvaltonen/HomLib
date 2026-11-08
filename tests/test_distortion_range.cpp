#include <Eigen/Dense>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <random>
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_kukelova_cvpr_2015.hpp"
#include "get_nakano_icpr_2025.hpp"
#include "get_valtonenornhag_icpr_2026.hpp"
#include "get_wadenback_3dv_2026.hpp"
#include "pose_estimator.h"
#include "posedata.hpp"

namespace {

HomLib::PoseData MakePose(double k1, double k2) {
    HomLib::PoseData pose;
    pose.homography = Eigen::Matrix3d::Identity();
    pose.distortion_parameter = k1;
    pose.distortion_parameter2 = k2;
    return pose;
}

// A solver that returns fixed solutions (some with NaN parameters), filtered by the range it
// gets, as the solvers do
class FixedSolver : public HomLib::PoseEstimator<FixedSolver> {
public:
    int solve(const std::vector<Eigen::Vector2d> &, const std::vector<Eigen::Vector2d> &,
              std::vector<HomLib::PoseData> *poses, const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        poses->push_back(MakePose(-0.3, -0.2));  // in the range
        poses->push_back(MakePose(nan, -0.2));
        poses->push_back(MakePose(-0.3, nan));
        poses->push_back(MakePose(nan, nan));
        poses->push_back(MakePose(-0.3, 0.0));   // one-sided
        poses->push_back(MakePose(-2.0, -0.2));  // out of the range
        HomLib::filter_distortion_range(poses, range);
        return poses->size();
    }
    int minimal_sample_size() const { return 5; }
};

bool HasK(const std::vector<HomLib::PoseData> &poses, double k, double tol) {
    for (const HomLib::PoseData &p : poses)
        if (std::abs(p.distortion_parameter - k) <= tol)
            return true;
    return false;
}

// A solver that records the range it gets
class RangeSolver : public HomLib::PoseEstimator<RangeSolver> {
public:
    mutable HomLib::DistortionRange last;
    int solve(const std::vector<Eigen::Vector2d> &, const std::vector<Eigen::Vector2d> &,
              std::vector<HomLib::PoseData> *poses, const HomLib::DistortionRange &range) const {
        last = range;
        poses->push_back(MakePose(-0.3, -0.2));
        return 1;
    }
    int minimal_sample_size() const { return 5; }
};

// Random two-sided problems in normalized coordinates (both images distorted), with a little
// noise, 6 matches
void MakeTwoSided(std::mt19937 &rng, std::vector<Eigen::Vector2d> &x, std::vector<Eigen::Vector2d> &y) {
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::normal_distribution<double> g(0.0, 1e-3);
    const double k1 = -0.1 * (u(rng) + 1), k2 = -0.1 * (u(rng) + 1);
    Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            H(i, j) += (i == 2 && j < 2 ? 0.05 : 0.2) * u(rng);
    auto distort = [](const Eigen::Vector2d &p, double k) {
        return Eigen::Vector2d(p * (2.0 / (1.0 + std::sqrt(1.0 - 4.0 * k * p.squaredNorm()))));
    };
    x.clear();
    y.clear();
    for (int i = 0; i < 6; i++) {
        const Eigen::Vector2d a(u(rng), u(rng));
        x.push_back(distort(a, k1) + Eigen::Vector2d(g(rng), g(rng)));
        y.push_back(distort((H * a.homogeneous()).hnormalized(), k2));
    }
}

bool Same(const HomLib::PoseData &a, const HomLib::PoseData &b) {
    // bit for bit (NaN entries included in the homography)
    return std::memcmp(a.homography.data(), b.homography.data(), 9 * sizeof(double)) == 0 &&
           a.distortion_parameter == b.distortion_parameter && a.distortion_parameter2 == b.distortion_parameter2;
}

}  // namespace

TEST_CASE("Distortion range: the check and its zero exemption") {
    const HomLib::DistortionRange off, on(-1.0, 0.1);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(off.contains(5.0, nan));
    REQUIRE(on.contains(-1.0, 0.1));
    REQUIRE(on.contains(0.0, -0.5));
    REQUIRE_FALSE(on.contains(-1.5, -0.5));
    REQUIRE_FALSE(on.contains(-0.5, 0.2));
    REQUIRE_FALSE(on.contains(nan, -0.5));
}

TEST_CASE("Distortion range: the estimators pass it to the solvers, in their normalized coordinates") {
    // The points don't change the range: it is in the solver's coordinates
    const std::vector<Eigen::Vector2d> x(5, Eigen::Vector2d(2, 2)), y = x;
    std::vector<HomLib::PoseData> poses;
    RangeSolver solver;
    solver.estimate(x, y, &poses);   // the default
    REQUIRE(solver.last.enabled);
    REQUIRE(solver.last.k_min == -1.0);
    REQUIRE(solver.last.k_max == 0.1);
    solver.k_range = HomLib::DistortionRange(-0.5, 0.05);
    solver.estimate(x, y, &poses);
    REQUIRE(solver.last.k_min == -0.5);
    REQUIRE(solver.last.k_max == 0.05);
    solver.k_range = HomLib::DistortionRange();
    solver.estimate(x, y, &poses);
    REQUIRE_FALSE(solver.last.enabled);
}

TEST_CASE("Distortion range in the solvers: the roots in the range, before their homography") {
    // The solvers that return all their roots: the same solutions as all of them, filtered
    const HomLib::DistortionRange range(-1.0, 0.1);
    namespace W = HomLib::Wadenback3DV2026;
    namespace K = HomLib::KukelovaCVPR2015;
    std::mt19937 rng(3);
    std::vector<Eigen::Vector2d> x, y;
    int checked = 0;
    for (int t = 0; t < 300; t++) {
        MakeTwoSided(rng, x, y);
        const std::vector<Eigen::Vector2d> x5(x.begin(), x.begin() + 5), y5(y.begin(), y.begin() + 5);
        std::vector<std::vector<HomLib::PoseData>> all = {
            W::get_double_sided(x5, y5, false), W::get_double_sided_equal(x5, y5, false), W::get_one_sided(x5, y5, false),
            K::get(x5, y5, false), K::get(x5, y5, true), K::get_6pt(x, y, false), K::get_6pt(x, y, true),
            HomLib::NakanoICPR2025::get(x5, y5, false)};
        const std::vector<std::vector<HomLib::PoseData>> in = {
            W::get_double_sided(x5, y5, false, range), W::get_double_sided_equal(x5, y5, false, range),
            W::get_one_sided(x5, y5, false, range), K::get(x5, y5, false, range), K::get(x5, y5, true, range),
            K::get_6pt(x, y, false, range), K::get_6pt(x, y, true, range),
            HomLib::NakanoICPR2025::get(x5, y5, false, range)};
        for (size_t s = 0; s < all.size(); s++) {
            HomLib::filter_distortion_range(&all[s], range);
            REQUIRE(all[s].size() == in[s].size());
            for (size_t i = 0; i < in[s].size(); i++)
                REQUIRE(Same(all[s][i], in[s][i]));
            checked += in[s].size();
        }
    }
    REQUIRE(checked > 300);
}

TEST_CASE("Distortion range in the solvers that keep one root: the best root in the range") {
    // With the check inside, these solvers choose among the roots in the range: a solution in the
    // range, the same as without the check when that one is in the range, and more often one
    const HomLib::DistortionRange range(-1.0, 0.1);
    namespace W = HomLib::Wadenback3DV2026;
    namespace F = HomLib::FitzgibbonCVPR2001;
    std::mt19937 rng(5);
    std::vector<Eigen::Vector2d> x, y;
    int gained = 0;
    for (int t = 0; t < 300; t++) {
        MakeTwoSided(rng, x, y);
        x.resize(5);
        y.resize(5);
        const std::vector<std::vector<HomLib::PoseData>> all = {
            W::get_double_sided_equal(x, y, true), W::get_one_sided(x, y, true), F::get(x, y), F::get_single(x, y),
            HomLib::NakanoICPR2025::get(x, y, true)};
        const std::vector<std::vector<HomLib::PoseData>> in = {
            W::get_double_sided_equal(x, y, true, range), W::get_one_sided(x, y, true, range), F::get(x, y, range),
            F::get_single(x, y, range), HomLib::NakanoICPR2025::get(x, y, true, range)};
        for (size_t s = 0; s < all.size(); s++) {
            REQUIRE(in[s].size() <= 1);
            for (const HomLib::PoseData &p : in[s])
                REQUIRE(range.contains(p.distortion_parameter, p.distortion_parameter2));
            std::vector<HomLib::PoseData> filtered = all[s];
            HomLib::filter_distortion_range(&filtered, range);
            if (!filtered.empty()) {
                REQUIRE(in[s].size() == 1);
                REQUIRE(Same(filtered[0], in[s][0]));
            }
            gained += filtered.empty() && !in[s].empty();
        }
    }
    REQUIRE(gained > 0);
}

TEST_CASE("Range check: NaN distortion parameters are dropped") {
    // Points with a normalization factor of 2: the distortion parameters of the input are those
    // of the solver divided by 4
    const std::vector<Eigen::Vector2d> x(5, Eigen::Vector2d(2, 2)), y = x;
    FixedSolver solver;
    std::vector<HomLib::PoseData> poses;
    // A range that excludes 0, so that the one-sided solution passes only as such
    solver.k_range = HomLib::DistortionRange(-1.0, -0.1);
    int n = solver.estimate(x, y, &poses);
    REQUIRE(n == 2);
    REQUIRE(poses.size() == 2);
    REQUIRE(poses[0].distortion_parameter == Catch::Approx(-0.3 / 4));
    REQUIRE(poses[0].distortion_parameter2 == Catch::Approx(-0.2 / 4));
    REQUIRE(poses[1].distortion_parameter == Catch::Approx(-0.3 / 4));
    REQUIRE(poses[1].distortion_parameter2 == 0.0);
}

TEST_CASE("Range check of the affine and orientation estimators") {
    // The problems of the Valtonen Ornhag ICPR 2026 tests, with the points multiplied by s: the
    // affine frames and orientations don't change, and the distortion parameter is divided by s^2
    const double s = 1000;
    std::vector<Eigen::Vector2d> xa = {
        Eigen::Vector2d(-0.304372059145248, -0.174125682172663),
        Eigen::Vector2d(-0.308082306270758,  0.100342577802908)
    };
    std::vector<Eigen::Vector2d> ya = {
        Eigen::Vector2d(0.596062339231747, 0.449851750693698),
        Eigen::Vector2d(0.284590553127076, 0.20125992440591)
    };
    const std::vector<Eigen::Matrix2d> A = {
        (Eigen::Matrix2d() << -0.742737063351745, -1.25720051435581,
                              1.29931162730836, -0.953395533218031).finished(),
        (Eigen::Matrix2d() << -0.763314941736403, -1.05311577674003,
                              1.13911319830585, -0.839438880618645).finished()
    };
    std::vector<Eigen::Vector2d> xo = {
        Eigen::Vector2d(0.611733207361777, -0.100815813483642),
        Eigen::Vector2d(0.506500688761625, 0.22654451870654),
        Eigen::Vector2d(0.618492431042353, 0.490634389234033),
        Eigen::Vector2d(0.446207436531052, -0.221543454776215)
    };
    std::vector<Eigen::Vector2d> yo = {
        Eigen::Vector2d(0.152378076324475, 0.471173894309057),
        Eigen::Vector2d(0.35249479911981, 0.604057206211103),
        Eigen::Vector2d(0.576429540272504, 0.546190419126011),
        Eigen::Vector2d(-0.0960825885072252, 0.57183541266364)
    };
    const std::vector<Eigen::Vector2d> ori = {
        Eigen::Vector2d(2.82616873462069, 2.11084467936878),
        Eigen::Vector2d(2.76670870212448, 1.5791745485192),
        Eigen::Vector2d(3.03912230291483, 1.40242855410085),
        Eigen::Vector2d(2.54398570702913, 2.11604504644583)
    };
    for (std::vector<Eigen::Vector2d> *v : {&xa, &ya, &xo, &yo})
        for (Eigen::Vector2d &q : *v)
            q *= s;
    const double ka = -0.15508855476035885 / (s * s), ko = -0.18221254202290402 / (s * s);
    const double tola = 1e-7 / (s * s);
    // The distortion parameters in the solvers' normalized coordinates
    auto norm2 = [](const std::vector<Eigen::Vector2d> &x) {
        double f0 = 0;
        for (const Eigen::Vector2d &q : x)
            f0 += q.norm();
        f0 /= x.size() * std::sqrt(2.0);
        return f0 * f0;
    };
    const double ka_norm = ka * norm2(xa), ko_norm = ko * norm2(xo);

    SECTION("affine") {
        HomLib::ValtonenOrnhagICPR2026::AffineSolverSingleSided solver;
        std::vector<HomLib::PoseData> poses;
        REQUIRE(solver.k_range.enabled);   // the default, [-1, 0.1]
        solver.estimate(xa, ya, A, &poses);
        REQUIRE(HasK(poses, ka, tola));

        // A range in the normalized coordinates, including or excluding k
        solver.k_range = HomLib::DistortionRange(2 * ka_norm, ka_norm / 2);
        solver.estimate(xa, ya, A, &poses);
        REQUIRE(HasK(poses, ka, tola));
        solver.k_range = HomLib::DistortionRange(ka_norm / 2, 0.1);
        REQUIRE(solver.estimate(xa, ya, A, &poses) == 0);
        REQUIRE(poses.empty());
    }

    SECTION("orientation") {
        HomLib::ValtonenOrnhagICPR2026::OrientationSolverSingleSided solver;
        std::vector<HomLib::PoseData> poses;
        REQUIRE(solver.k_range.enabled);   // the default, [-1, 0.1]
        solver.estimate(xo, yo, ori, &poses);
        REQUIRE(HasK(poses, ko, tola));

        // A range in the normalized coordinates, including or excluding k
        solver.k_range = HomLib::DistortionRange(2 * ko_norm, ko_norm / 2);
        solver.estimate(xo, yo, ori, &poses);
        REQUIRE(HasK(poses, ko, tola));
        solver.k_range = HomLib::DistortionRange(ko_norm / 2, 0.1);
        REQUIRE(solver.estimate(xo, yo, ori, &poses) == 0);
        REQUIRE(poses.empty());
    }
}
