#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "posedata.hpp"
#include "problem_instance.hpp"
#include "refinement.hpp"

namespace {

Eigen::Vector2d Distort(const Eigen::Vector2d &u, double k)
{
    return u * (2.0 / (1.0 + std::sqrt(1.0 - 4.0 * k * u.squaredNorm())));
}

// The distortion cases of the refinements (HomLib::DistortionCase): NO_DISTORTION
// (refinement_no_dist), ONE_SIDED_LEFT, x2 only (refinement_onesided, k = distortion_parameter2),
// ONE_SIDED_RIGHT, x1 only (refinement_onesided_right, k = distortion_parameter), TWO_SIDED_EQUAL
// (refinement_twosided_equal) and TWO_SIDED (refinement_twosided)
using HomLib::DistortionCase;

// The cost minimized by the refinements, at q = (H entries in column-major order, h22 = 1; k1; k2),
// with k1 = 0 or k2 = 0 where the model has no distortion: for the two-sided models, the transfer
// error in the second image, sum |Distort(pi(H (x1, 1 + k1 |x1|^2)), k2) - x2|^2; for the others,
// sum |pi(H (x1, 1 + k1 |x1|^2)) - x2 / (1 + k2 |x2|^2)|^2
double Cost(DistortionCase model, const std::vector<Eigen::Vector2d> &x1, const std::vector<Eigen::Vector2d> &x2,
            const double *q)
{
    Eigen::Matrix3d H;
    H << q[0], q[3], q[6], q[1], q[4], q[7], q[2], q[5], 1.0;
    const bool d1 = model == DistortionCase::ONE_SIDED_RIGHT || model == DistortionCase::TWO_SIDED_EQUAL || model == DistortionCase::TWO_SIDED;
    const bool d2 = model == DistortionCase::ONE_SIDED_LEFT || model == DistortionCase::TWO_SIDED_EQUAL || model == DistortionCase::TWO_SIDED;
    const double k1 = d1 ? q[8] : 0.0, k2 = model == DistortionCase::TWO_SIDED_EQUAL ? q[8] : d2 ? q[9] : 0.0;
    double c = 0;
    for (size_t i = 0; i < x1.size(); i++) {
        const Eigen::Vector2d z = (H * Eigen::Vector3d(x1[i](0), x1[i](1), 1.0 + k1 * x1[i].squaredNorm())).hnormalized();
        if (model == DistortionCase::TWO_SIDED_EQUAL || model == DistortionCase::TWO_SIDED)
            c += (Distort(z, k2) - x2[i]).squaredNorm();
        else
            c += (z - x2[i] / (1.0 + k2 * x2[i].squaredNorm())).squaredNorm();
    }
    return c;
}

// The largest partial derivative of the cost, by central differences, relative to the cost
double RelativeGradient(DistortionCase model, const std::vector<Eigen::Vector2d> &x1, const std::vector<Eigen::Vector2d> &x2,
                        const HomLib::PoseData &p)
{
    const Eigen::Matrix3d H = p.homography / p.homography(2, 2);
    double q[10] = {H(0, 0), H(1, 0), H(2, 0), H(0, 1), H(1, 1), H(2, 1), H(0, 2), H(1, 2), p.distortion_parameter,
                    p.distortion_parameter2};
    std::vector<int> params = {0, 1, 2, 3, 4, 5, 6, 7};
    if (model == DistortionCase::ONE_SIDED_RIGHT || model == DistortionCase::TWO_SIDED_EQUAL || model == DistortionCase::TWO_SIDED)
        params.push_back(8);
    if (model == DistortionCase::ONE_SIDED_LEFT || model == DistortionCase::TWO_SIDED)
        params.push_back(9);
    double worst = 0;
    for (int j : params) {
        const double h = 1e-7, qj = q[j];
        q[j] = qj + h;
        const double cp = Cost(model, x1, x2, q);
        q[j] = qj - h;
        const double cm = Cost(model, x1, x2, q);
        q[j] = qj;
        worst = std::max(worst, std::abs(cp - cm) / (2 * h));
    }
    return worst / Cost(model, x1, x2, q);
}

}  // namespace

TEST_CASE("Refinements: the gradient of the cost vanishes at the result, for every distortion case") {
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::normal_distribution<double> g(0.0, 1e-3);
    for (DistortionCase model : {DistortionCase::NO_DISTORTION, DistortionCase::ONE_SIDED_LEFT, DistortionCase::ONE_SIDED_RIGHT, DistortionCase::TWO_SIDED_EQUAL, DistortionCase::TWO_SIDED}) {
        const bool d1 = model == DistortionCase::ONE_SIDED_RIGHT || model == DistortionCase::TWO_SIDED_EQUAL || model == DistortionCase::TWO_SIDED;
        const bool d2 = model == DistortionCase::ONE_SIDED_LEFT || model == DistortionCase::TWO_SIDED_EQUAL || model == DistortionCase::TWO_SIDED;
        double worst = 0;
        for (int t = 0; t < 20; t++) {
            double k1 = -0.2 + 0.05 * u(rng), k2 = model == DistortionCase::TWO_SIDED ? -0.1 + 0.05 * u(rng) : k1;
            if (!d1)
                k1 = 0.0;
            if (!d2)
                k2 = 0.0;
            Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    H(i, j) += (i == 2 && j < 2 ? 0.05 : 0.2) * u(rng);
            std::vector<Eigen::Vector2d> x1, x2;
            for (int i = 0; i < 100; i++) {
                const Eigen::Vector2d a(u(rng), u(rng));
                const Eigen::Vector2d z = (H * Eigen::Vector3d(a(0), a(1), 1.0 + k1 * a.squaredNorm())).hnormalized();
                x1.push_back(a);
                x2.push_back((d2 ? Distort(z, k2) : z) + Eigen::Vector2d(g(rng), g(rng)));
            }
            HomLib::PoseData p;
            p.homography = H / H(2, 2);
            p.distortion_parameter = k1;
            p.distortion_parameter2 = k2;
            if (model == DistortionCase::NO_DISTORTION)
                HomLib::refinement_no_dist(x1, x2, p);
            else if (model == DistortionCase::ONE_SIDED_LEFT)
                HomLib::refinement_onesided(x1, x2, p);
            else if (model == DistortionCase::ONE_SIDED_RIGHT)
                HomLib::refinement_onesided_right(x1, x2, p);
            else if (model == DistortionCase::TWO_SIDED_EQUAL)
                HomLib::refinement_twosided_equal(x1, x2, p);
            else
                HomLib::refinement_twosided(x1, x2, p);
            worst = std::max(worst, RelativeGradient(model, x1, x2, p));
        }
        CAPTURE(model, worst);
        CHECK(worst < 1e-4);
    }
}

TEST_CASE("Refinements: a model of an identity with distortion still fits its points") {
    // With H = I and k1 = k2 = t, the model maps every point of the first image to itself (where
    // 4 t |x|^2 < 1): t isn't observable. The two-sided refinements minimized the error in the
    // undistorted plane of the second image, which shrinks as t grows: from the ground truth, they
    // went to t ~ 100 and errors of a third of the points' spread (HPatches' illumination pairs,
    // with every two-sided solver, on GC-RANSAC and MAGSAC). They now minimize the error in the
    // second image, so the refined model fits the points (t can still move). The one-sided
    // refinement keeps its cost, which needs H to collapse as well to tend to 0.
    std::mt19937 rng(2);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::normal_distribution<double> g(0.0, 3e-3);
    std::vector<Eigen::Vector2d> x1, x2;
    for (int i = 0; i < 300; i++) {
        x1.emplace_back(u(rng), u(rng));
        x2.push_back(x1.back() + Eigen::Vector2d(g(rng), g(rng)));
    }
    for (DistortionCase model : {DistortionCase::ONE_SIDED_LEFT, DistortionCase::TWO_SIDED_EQUAL, DistortionCase::TWO_SIDED}) {
        HomLib::PoseData p;
        p.homography = Eigen::Matrix3d::Identity();
        p.distortion_parameter = p.distortion_parameter2 = 0.0;
        if (model == DistortionCase::ONE_SIDED_LEFT)
            HomLib::refinement_onesided(x1, x2, p);
        else if (model == DistortionCase::TWO_SIDED_EQUAL)
            HomLib::refinement_twosided_equal(x1, x2, p);
        else
            HomLib::refinement_twosided(x1, x2, p);
        const double k1 = model == DistortionCase::ONE_SIDED_LEFT ? 0.0 : p.distortion_parameter;
        const double k2 = p.distortion_parameter2;
        std::vector<double> errors;
        for (size_t i = 0; i < x1.size(); i++) {
            const Eigen::Vector2d z = (p.homography * Eigen::Vector3d(x1[i](0), x1[i](1), 1.0 + k1 * x1[i].squaredNorm())).hnormalized();
            errors.push_back(4 * k2 * z.squaredNorm() < 1 ? (Distort(z, k2) - x2[i]).norm() : 1e9);
        }
        std::nth_element(errors.begin(), errors.begin() + errors.size() / 2, errors.end());
        const double median = errors[errors.size() / 2];
        CAPTURE(model, k1, k2, median);
        CHECK(median < 5e-3);
    }
}
