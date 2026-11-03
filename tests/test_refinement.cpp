#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "posedata.hpp"
#include "refinement.hpp"

namespace {

// Point of the division model with parameter k whose undistorted point is u
Eigen::Vector2d Distort(const Eigen::Vector2d &u, double k)
{
    return u * (2.0 / (1.0 + std::sqrt(1.0 - 4.0 * k * u.squaredNorm())));
}

enum Model { ONE_SIDED_RIGHT, TWO_SIDED_EQUAL, TWO_SIDED };

// Sum of the squared residuals of the refinement functions: z = pi(H (x1, 1 + k1 |x1|^2))
// against x2 (one-sided right) or x2 / (1 + k2 |x2|^2) (two-sided; k2 = k1 if equal), for the
// parameters q = (H in column-major order without H(2, 2) = 1, k1[, k2])
double Cost(Model model, const std::vector<Eigen::Vector2d> &x1, const std::vector<Eigen::Vector2d> &x2,
            const double *q)
{
    Eigen::Matrix3d H;
    H << q[0], q[3], q[6], q[1], q[4], q[7], q[2], q[5], 1.0;
    const double k1 = q[8], k2 = model == TWO_SIDED ? q[9] : q[8];
    double c = 0;
    for (size_t i = 0; i < x1.size(); i++) {
        const Eigen::Vector2d z = (H * Eigen::Vector3d(x1[i](0), x1[i](1), 1.0 + k1 * x1[i].squaredNorm())).hnormalized();
        const Eigen::Vector2d y = model == ONE_SIDED_RIGHT ? x2[i] : Eigen::Vector2d(x2[i] / (1.0 + k2 * x2[i].squaredNorm()));
        c += (z - y).squaredNorm();
    }
    return c;
}

// The largest component of the cost's gradient (central differences), relative to the cost
double RelativeGradient(Model model, const std::vector<Eigen::Vector2d> &x1, const std::vector<Eigen::Vector2d> &x2,
                        const HomLib::PoseData &p)
{
    const Eigen::Matrix3d H = p.homography / p.homography(2, 2);
    double q[10] = {H(0, 0), H(1, 0), H(2, 0), H(0, 1), H(1, 1), H(2, 1), H(0, 2), H(1, 2), p.distortion_parameter,
                    p.distortion_parameter2};
    const int n = model == TWO_SIDED ? 10 : 9;
    double worst = 0;
    for (int j = 0; j < n; j++) {
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

// At the result of the refinement, the gradient of its cost vanishes (a least-squares
// optimum): with the Jacobian of the residuals right, including the derivatives with respect
// to H(0, 2) and H(1, 2), which multiply 1 + k1 |x1|^2
TEST_CASE("Refinement with the distortion of x1: the gradient of the cost vanishes at the result") {
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::normal_distribution<double> g(0.0, 1e-3);
    for (Model model : {ONE_SIDED_RIGHT, TWO_SIDED_EQUAL, TWO_SIDED}) {
        double worst = 0;
        for (int t = 0; t < 20; t++) {
            const double k1 = -0.2 + 0.05 * u(rng), k2 = model == TWO_SIDED ? -0.1 + 0.05 * u(rng) : k1;
            Eigen::Matrix3d H = Eigen::Matrix3d::Identity();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    H(i, j) += (i == 2 && j < 2 ? 0.05 : 0.2) * u(rng);
            std::vector<Eigen::Vector2d> x1, x2;
            for (int i = 0; i < 100; i++) {
                const Eigen::Vector2d a(u(rng), u(rng));
                const Eigen::Vector2d z = (H * Eigen::Vector3d(a(0), a(1), 1.0 + k1 * a.squaredNorm())).hnormalized();
                x1.push_back(a);
                x2.push_back((model == ONE_SIDED_RIGHT ? z : Distort(z, k2)) + Eigen::Vector2d(g(rng), g(rng)));
            }
            HomLib::PoseData p;
            p.homography = H / H(2, 2);
            p.distortion_parameter = k1;
            p.distortion_parameter2 = model == TWO_SIDED ? k2 : 0.0;
            if (model == ONE_SIDED_RIGHT)
                HomLib::refinement_onesided_right(x1, x2, p);
            else if (model == TWO_SIDED_EQUAL)
                HomLib::refinement_twosided_equal(x1, x2, p);
            else
                HomLib::refinement_twosided(x1, x2, p);
            worst = std::max(worst, RelativeGradient(model, x1, x2, p));
        }
        CAPTURE(model, worst);
        REQUIRE(worst < 1e-4);
    }
}
