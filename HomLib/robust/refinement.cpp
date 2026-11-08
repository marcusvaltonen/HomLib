#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include "posedata.hpp"

static const double SMALL_NUMBER = 1e-8;
static const double TOL_CONVERGENCE = 1e-10;
static const double INITIAL_LM_DAMP = 1e-4;
static const int MAX_ITER = 10;

namespace {
// Levenberg-Marquardt on H (H(2,2) = 1) and the distortion parameters of p, the loop that the
// refinements below share: fill(H, J, res) computes the residuals and their Jacobian (in the
// order h11, h21, h31, h12, h22, h32, h13, h23, then the distortion parameters) at H and p's
// distortion parameters, and update_dist(dx) adds the steps of the distortion parameters to p.
template <class Fill, class UpdateDist>
void refine(size_t n_pts, int n_params, HomLib::PoseData &p, const Fill &fill, const UpdateDist &update_dist) {
    double lm_damp = INITIAL_LM_DAMP;
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2 * n_pts, n_params);
    Eigen::VectorXd res = Eigen::VectorXd::Zero(2 * n_pts);
    for (int iter = 0; iter < MAX_ITER; iter++) {
        Eigen::Matrix3d H = p.homography;
        H /= H(2, 2);
        fill(H, J, res);
        if (res.norm() < TOL_CONVERGENCE)
            break;
        Eigen::MatrixXd Hess = J.transpose() * J;
        Hess.diagonal().array() += lm_damp;   // LM dampening
        const Eigen::VectorXd g = -J.transpose() * res;
        if (g.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
            break;
        const Eigen::VectorXd dx = Hess.ldlt().solve(g);
        if (!dx.allFinite())   // a degenerate model or sample: keep the last finite model
            break;
        Eigen::Map<Eigen::Matrix<double, 8, 1>>(H.data()) += dx.head(8);
        p.homography = H;
        update_dist(dx);
        if (dx.array().abs().maxCoeff() < SMALL_NUMBER)
            break;
        lm_damp = std::max(1e-8, lm_damp / 10.0);
    }
}
}  // namespace

namespace HomLib {
// Non-linear refinement of transfer error |x2 - pi(H*x1)|^2, parameterized by fixing H(2,2) = 1
void refinement_no_dist(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    refine(x1.size(), 8, p, [&](const Eigen::Matrix3d &H, Eigen::MatrixXd &J, Eigen::VectorXd &res) {
        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);

        for (size_t i = 0; i < x1.size(); i++) {
            const double x1_0 = x1[i](0);
            const double x1_1 = x1[i](1);
            const double x2_0 = x2[i](0);
            const double x2_1 = x2[i](1);

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;

            // Residual
            res(2 * i + 0) = z0 - x2_0;
            res(2 * i + 1) = z1 - x2_1;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = 1.0 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;

            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = 1.0 * inv_Hx1_2;

        }

    }, [](const Eigen::VectorXd &) {});
}
// Non-linear refinement of transfer error |x2(k) - pi(H*x1)|^2, parameterized by fixing H(2,2) = 1
void refinement_onesided(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    refine(x1.size(), 9, p, [&](const Eigen::Matrix3d &H, Eigen::MatrixXd &J, Eigen::VectorXd &res) {
        const double k = p.distortion_parameter2;

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);

        for (size_t i = 0; i < x1.size(); i++) {
            const double x1_0 = x1[i](0);
            const double x1_1 = x1[i](1);
            const double x2_0 = x2[i](0);
            const double x2_1 = x2[i](1);

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;

            // Residual
            const double r22 = (x2_0 * x2_0 + x2_1 * x2_1);
            const double dist_factor = 1.0 / (1.0 + k * r22);
            res(2 * i + 0) = z0 - x2_0 * dist_factor;
            res(2 * i + 1) = z1 - x2_1 * dist_factor;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = 1.0 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;
            J(2 * i + 0, 8) = x2_0 * r22 * dist_factor * dist_factor;  // Dist coeff.

            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = 1.0 * inv_Hx1_2;
            J(2 * i + 1, 8) = x2_1 * r22 * dist_factor * dist_factor;  // Dist coeff.

        }

    }, [&p](const Eigen::VectorXd &dx) {
        p.distortion_parameter2 += dx(8);
    });
}

// Non-linear refinement of transfer error |x2 - pi(H*x1(k))|^2, parameterized by fixing H(2,2) = 1
void refinement_onesided_right(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    refine(x1.size(), 9, p, [&](const Eigen::Matrix3d &H, Eigen::MatrixXd &J, Eigen::VectorXd &res) {
        const double k = p.distortion_parameter;

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);

        for (size_t i = 0; i < x1.size(); i++) {          
            double x1_0 = x1[i](0);
            double x1_1 = x1[i](1);
            double x2_0 = x2[i](0);
            double x2_1 = x2[i](1);

            const double r12 = (x1_0 * x1_0 + x1_1 * x1_1);
            const double inv_dist_factor1 = 1.0 + k * r12;

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2 * inv_dist_factor1;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2 * inv_dist_factor1;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2 * inv_dist_factor1);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;

            // Residual
            res(2 * i + 0) = z0 - x2_0;
            res(2 * i + 1) = z1 - x2_1;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = inv_dist_factor1 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;
            J(2 * i + 0, 8) = H0_2 * r12 * inv_Hx1_2 - z0 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff. k

            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = inv_dist_factor1 * inv_Hx1_2;
            J(2 * i + 1, 8) = H1_2 * r12 * inv_Hx1_2 - z1 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff. k
        }

    }, [&p](const Eigen::VectorXd &dx) {
        p.distortion_parameter += dx(8);
    });
}
// Non-linear refinement of transfer error |x2(k) - pi(H*x1(k))|^2, parameterized by fixing H(2,2) = 1
void refinement_twosided_equal(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    refine(x1.size(), 9, p, [&](const Eigen::Matrix3d &H, Eigen::MatrixXd &J, Eigen::VectorXd &res) {
        const double k = p.distortion_parameter;

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);

        for (size_t i = 0; i < x1.size(); i++) {          
            double x1_0 = x1[i](0);
            double x1_1 = x1[i](1);
            double x2_0 = x2[i](0);
            double x2_1 = x2[i](1);

            const double r12 = (x1_0 * x1_0 + x1_1 * x1_1);
            const double inv_dist_factor1 = 1.0 + k * r12;

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2 * inv_dist_factor1;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2 * inv_dist_factor1;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2 * inv_dist_factor1);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;

            // Residual
            const double r22 = (x2_0 * x2_0 + x2_1 * x2_1);
            const double dist_factor2 = 1.0 / (1.0 + k * r22);
            res(2 * i + 0) = z0 - x2_0 * dist_factor2;
            res(2 * i + 1) = z1 - x2_1 * dist_factor2;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = inv_dist_factor1 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;
            J(2 * i + 0, 8) = x2_0 * r22 * dist_factor2 * dist_factor2 + H0_2 * r12 * inv_Hx1_2 - z0 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff.

            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = inv_dist_factor1 * inv_Hx1_2;
            J(2 * i + 1, 8) = x2_1 * r22 * dist_factor2 * dist_factor2 + H1_2 * r12 * inv_Hx1_2 - z1 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff.
        }

    }, [&p](const Eigen::VectorXd &dx) {
        p.distortion_parameter += dx(8);
        p.distortion_parameter2 = p.distortion_parameter;
    });
}
// Non-linear refinement of transfer error |x2(k2) - pi(H*x1(k1))|^2, parameterized by fixing H(2,2) = 1
void refinement_twosided(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    refine(x1.size(), 10, p, [&](const Eigen::Matrix3d &H, Eigen::MatrixXd &J, Eigen::VectorXd &res) {
        const double k1 = p.distortion_parameter;
        const double k2 = p.distortion_parameter2;

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);

        for (size_t i = 0; i < x1.size(); i++) {          
            double x1_0 = x1[i](0);
            double x1_1 = x1[i](1);
            double x2_0 = x2[i](0);
            double x2_1 = x2[i](1);

            const double r12 = (x1_0 * x1_0 + x1_1 * x1_1);
            const double inv_dist_factor1 = 1.0 + k1 * r12;

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2 * inv_dist_factor1;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2 * inv_dist_factor1;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2 * inv_dist_factor1);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;

            // Residual
            const double r22 = (x2_0 * x2_0 + x2_1 * x2_1);
            const double dist_factor2 = 1.0 / (1.0 + k2 * r22);
            res(2 * i + 0) = z0 - x2_0 * dist_factor2;
            res(2 * i + 1) = z1 - x2_1 * dist_factor2;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = inv_dist_factor1 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;
            J(2 * i + 0, 8) = H0_2 * r12 * inv_Hx1_2 - z0 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff. k1
            J(2 * i + 0, 9) = x2_0 * r22 * dist_factor2 * dist_factor2;  // Dist coeff. k2

            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = inv_dist_factor1 * inv_Hx1_2;
            J(2 * i + 1, 8) = H1_2 * r12 * inv_Hx1_2 - z1 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff. k1
            J(2 * i + 1, 9) = x2_1 * r22 * dist_factor2 * dist_factor2;  // Dist coeff. k2
        }

    }, [&p](const Eigen::VectorXd &dx) {
        p.distortion_parameter += dx(8);
        p.distortion_parameter2 += dx(9);
    });
}
}

