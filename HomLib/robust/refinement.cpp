#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <vector>
#include "posedata.hpp"

static const double SMALL_NUMBER = 1e-8;
static const double TOL_CONVERGENCE = 1e-10;
static const double TOL_CONVERGENCE_SQ = TOL_CONVERGENCE * TOL_CONVERGENCE;
static const double INITIAL_LM_DAMP = 1e-4;
static const int MAX_ITER = 10;

namespace {

// Levenberg-Marquardt on H (H(2,2) = 1) and the distortion parameters of p, shared by the
// refinements below. accumulate(H, JtJ, g, cost) adds every point's contribution to JtJ (upper
// triangle), g = -J^T r and the cost, with the parameters in the order h11, h21, h31, h12, h22,
// h32, h13, h23, then the distortion parameters, and returns false if a point has no residual at
// p (see accumulate_distorted); update_dist(dx) adds the distortion steps to p. A step is kept
// only if it lowers the cost (except a last step below SMALL_NUMBER); otherwise the damping grows
// tenfold and the step is recomputed. At most MAX_ITER steps are kept, and MAX_ITER rejected.
template <int N, class Accumulate, class UpdateDist>
void refine(HomLib::PoseData &p, const Accumulate &accumulate, const UpdateDist &update_dist) {
    double lm_damp = INITIAL_LM_DAMP;
    HomLib::PoseData best = p;
    Eigen::Matrix<double, N, N> JtJ;
    Eigen::Matrix<double, N, 1> g;
    double best_cost = 0.0;
    for (int kept = 0, rejected = 0; kept <= MAX_ITER && rejected <= MAX_ITER;) {
        Eigen::Matrix3d H = p.homography;
        H /= H(2, 2);
        Eigen::Matrix<double, N, N> JtJ_p = Eigen::Matrix<double, N, N>::Zero();
        Eigen::Matrix<double, N, 1> g_p = Eigen::Matrix<double, N, 1>::Zero();
        double cost = 0.0;
        const bool ok = accumulate(H, JtJ_p, g_p, cost);
        if (kept == 0 && rejected == 0 && !(ok && std::isfinite(cost)))
            return;   // no residual at the initial model: keep it
        if (kept == 0 || (ok && cost < best_cost)) {
            if (kept > 0)
                lm_damp = std::max(1e-8, lm_damp / 10.0);
            best = p;
            best.homography = H;
            JtJ = JtJ_p;
            g = g_p;
            best_cost = cost;
            kept++;
        } else {
            lm_damp *= 10.0;
            rejected++;
        }
        if (kept > MAX_ITER || best_cost < TOL_CONVERGENCE_SQ || g.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
            break;

        Eigen::Matrix<double, N, N> A = JtJ;
        A.diagonal().array() += lm_damp;
        const Eigen::Matrix<double, N, 1> dx = A.template selfadjointView<Eigen::Upper>().ldlt().solve(g);
        if (!dx.allFinite())   // a degenerate model or sample: keep the last finite model
            break;
        p = best;
        Eigen::Map<Eigen::Matrix<double, 8, 1>>(p.homography.data()) += dx.head(8);
        update_dist(dx);
        if (dx.array().abs().maxCoeff() < SMALL_NUMBER)
            return;   // converged: the last step, too small to be checked
    }
    p = best;
}

// Accumulate the H-parameter block (cols 0-7) of JtJ (upper triangle) and g = -J^T r, for one
// point's two Jacobian rows:
//   even row: [a, 0, -az0, b, 0, -bz0, c, 0]
//   odd row:  [0, a, -az1, 0, b, -bz1, 0, c]
// (az0, bz0, az1, bz1 are passed positive; the Jacobian entries are their negatives.)
//
// The even/odd Jacobian rows have complementary sparsity: cols {0,3,6} are
// zero in odd rows, cols {1,4,7} are zero in even rows.  This gives 9
// structural zeros in the 8x8 block and makes the {0,3,6} diagonal block
// identical to the {1,4,7} block (AA = BB), so 6 products serve both.
template <int N>
inline void accumulate_h_block(
    Eigen::Matrix<double, N, N> &JtJ, Eigen::Matrix<double, N, 1> &g, double &cost,
    double a, double b, double c,
    double az0, double bz0, double az1, double bz1,
    double re, double ro)
{
    cost += re * re + ro * ro;

    const double aa = a * a, ab = a * b, ac = a * c;
    const double bb = b * b, bc = b * c, cc = c * c;
    JtJ(0, 0) += aa;  JtJ(1, 1) += aa;
    JtJ(0, 3) += ab;  JtJ(1, 4) += ab;
    JtJ(0, 6) += ac;  JtJ(1, 7) += ac;
    JtJ(3, 3) += bb;  JtJ(4, 4) += bb;
    JtJ(3, 6) += bc;  JtJ(4, 7) += bc;
    JtJ(6, 6) += cc;  JtJ(7, 7) += cc;

    JtJ(0, 2) -= a * az0;
    JtJ(0, 5) -= a * bz0;
    JtJ(2, 3) -= az0 * b;
    JtJ(3, 5) -= b * bz0;
    JtJ(2, 6) -= az0 * c;
    JtJ(5, 6) -= bz0 * c;

    JtJ(1, 2) -= a * az1;
    JtJ(1, 5) -= a * bz1;
    JtJ(2, 4) -= az1 * b;
    JtJ(4, 5) -= b * bz1;
    JtJ(2, 7) -= az1 * c;
    JtJ(5, 7) -= bz1 * c;

    JtJ(2, 2) += az0 * az0 + az1 * az1;
    JtJ(2, 5) += az0 * bz0 + az1 * bz1;
    JtJ(5, 5) += bz0 * bz0 + bz1 * bz1;

    g(0) -= a * re;
    g(1) -= a * ro;
    g(2) += az0 * re + az1 * ro;
    g(3) -= b * re;
    g(4) -= b * ro;
    g(5) += bz0 * re + bz1 * ro;
    g(6) -= c * re;
    g(7) -= c * ro;
}

// Accumulate a single distortion column D into JtJ (upper triangle) and g.
// dk_e / dk_o are the even / odd Jacobian entries for column D.
template <int N, int D>
inline void accumulate_dist_col(
    Eigen::Matrix<double, N, N> &JtJ, Eigen::Matrix<double, N, 1> &g,
    double a, double b, double c,
    double az0, double bz0, double az1, double bz1,
    double dk_e, double dk_o, double re, double ro)
{
    JtJ(0, D) += a * dk_e;
    JtJ(3, D) += b * dk_e;
    JtJ(6, D) += c * dk_e;
    JtJ(1, D) += a * dk_o;
    JtJ(4, D) += b * dk_o;
    JtJ(7, D) += c * dk_o;
    JtJ(2, D) -= az0 * dk_e + az1 * dk_o;
    JtJ(5, D) -= bz0 * dk_e + bz1 * dk_o;
    JtJ(D, D) += dk_e * dk_e + dk_o * dk_o;
    g(D) -= dk_e * re + dk_o * ro;
}

// The residual in the second image, where its points are distorted with k2: r = d(z) - x2, where
// z = pi(H x1(k1)) is the first image's point mapped to the second image's undistorted plane, and
// d(z) = 2 z / (1 + sqrt(1 - 4 k2 |z|^2)) distorts it (the inverse of x2 / (1 + k2 |x2|^2), the
// root that tends to z as k2 -> 0, as in try_radialdistort). The residual is the backends' and
// the HPatches benchmark's, in pixels of the second image. (The residual in the undistorted plane,
// z - x2 / (1 + k2 |x2|^2), tends to 0 for any matches as k1 = k2 grows with h13 = h23 = 0: the
// refinement drifted there from the true model of an identity.) jz0 and jz1: the derivatives of
// z0 and z1 with respect to the N parameters; k2_col: the column of k2 (which may also be k1's).
// Adds the point's contribution as accumulate_h_block does; false, adding nothing, if z can't be
// distorted (4 k2 |z|^2 >= 1).
template <int N>
inline bool accumulate_distorted(
    Eigen::Matrix<double, N, N> &JtJ, Eigen::Matrix<double, N, 1> &g, double &cost,
    double z0, double z1, double k2, double p0, double p1,
    const Eigen::Matrix<double, N, 1> &jz0, const Eigen::Matrix<double, N, 1> &jz1, int k2_col)
{
    const double s = z0 * z0 + z1 * z1;
    const double q2 = 1.0 - 4.0 * k2 * s;
    if (!(q2 > 0.0) || !std::isfinite(s))   // also false for a NaN
        return false;
    const double q = std::sqrt(q2);
    const double f = 2.0 / (1.0 + q);
    // df/ds = k2 h and df/dk2 = s h
    const double h = 4.0 / (q * (1.0 + q) * (1.0 + q));
    const double re = f * z0 - p0, ro = f * z1 - p1;
    cost += re * re + ro * ro;
    // d(d(z))/dz = f I + 2 k2 h z z^T
    const double m00 = f + 2.0 * k2 * h * z0 * z0, m01 = 2.0 * k2 * h * z0 * z1, m11 = f + 2.0 * k2 * h * z1 * z1;
    Eigen::Matrix<double, N, 1> je = m00 * jz0 + m01 * jz1, jo = m01 * jz0 + m11 * jz1;
    je(k2_col) += s * h * z0;
    jo(k2_col) += s * h * z1;
    JtJ.template selfadjointView<Eigen::Upper>().rankUpdate(je);
    JtJ.template selfadjointView<Eigen::Upper>().rankUpdate(jo);
    g -= je * re + jo * ro;
    return true;
}

// The derivatives of z = pi(H (u, v, d1)) with respect to the eight entries of H, in the order of
// refine, given a = u w, b = v w, c = d1 w and z, where w = 1 / (h31 u + h32 v + h33 d1)
template <int N>
inline void h_derivatives(Eigen::Matrix<double, N, 1> &jz0, Eigen::Matrix<double, N, 1> &jz1,
                          double a, double b, double c, double z0, double z1)
{
    jz0.setZero();
    jz1.setZero();
    jz0(0) = a; jz0(2) = -a * z0; jz0(3) = b; jz0(5) = -b * z0; jz0(6) = c;
    jz1(1) = a; jz1(2) = -a * z1; jz1(4) = b; jz1(5) = -b * z1; jz1(7) = c;
}

// Non-linear refinement of transfer error |d(pi(H*x1(k1)), k2) - x2|^2, where x1(k1) undistorts
// x1 with k1 and d distorts with k2 (see accumulate_distorted). Equal: k1 = k2.
template <bool Equal>
void refinement_twosided_impl(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p)
{
    constexpr int N = Equal ? 9 : 10;
    refine<N>(p,
        [&](const Eigen::Matrix3d &H, auto &JtJ, auto &g, double &cost) {
            const double k1 = p.distortion_parameter;
            const double k2 = Equal ? k1 : p.distortion_parameter2;
            Eigen::Matrix<double, N, 1> jz0, jz1;
            for (size_t i = 0; i < x1.size(); i++) {
                const double u = x1[i](0), v = x1[i](1);
                const double r12 = u * u + v * v;
                const double d1 = 1.0 + k1 * r12;
                const double w = 1.0 / (H(2, 0) * u + H(2, 1) * v + H(2, 2) * d1);
                const double z0 = (H(0, 0) * u + H(0, 1) * v + H(0, 2) * d1) * w;
                const double z1 = (H(1, 0) * u + H(1, 1) * v + H(1, 2) * d1) * w;
                h_derivatives<N>(jz0, jz1, u * w, v * w, d1 * w, z0, z1);
                jz0(8) = r12 * w * (H(0, 2) - z0 * H(2, 2));
                jz1(8) = r12 * w * (H(1, 2) - z1 * H(2, 2));
                if (!accumulate_distorted<N>(JtJ, g, cost, z0, z1, k2, x2[i](0), x2[i](1), jz0, jz1, N - 1))
                    return false;
            }
            return true;
        },
        [&p](const auto &dx) {
            p.distortion_parameter += dx(8);
            p.distortion_parameter2 = Equal ? p.distortion_parameter : p.distortion_parameter2 + dx(N - 1);
        });
}

}  // namespace

namespace HomLib {

// Non-linear refinement of transfer error |x2 - pi(H*x1)|^2, parameterized by fixing H(2,2) = 1
void refinement_no_dist(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p)
{
    refine<8>(p,
        [&](const Eigen::Matrix3d &H, auto &JtJ, auto &g, double &cost) {
            const double H00 = H(0, 0), H01 = H(0, 1), H02 = H(0, 2);
            const double H10 = H(1, 0), H11 = H(1, 1), H12 = H(1, 2);
            const double H20 = H(2, 0), H21 = H(2, 1);

            for (size_t i = 0; i < x1.size(); i++) {
                const double u = x1[i](0), v = x1[i](1);
                const double p0 = x2[i](0), p1 = x2[i](1);

                const double w = 1.0 / (H20 * u + H21 * v + H(2, 2));
                const double z0 = (H00 * u + H01 * v + H02) * w;
                const double z1 = (H10 * u + H11 * v + H12) * w;

                accumulate_h_block(JtJ, g, cost,
                    u * w, v * w, w,
                    u * w * z0, v * w * z0, u * w * z1, v * w * z1,
                    z0 - p0, z1 - p1);
            }
            return true;
        },
        [](const auto &) {});
}

// Non-linear refinement of transfer error |x2(k) - pi(H*x1)|^2
void refinement_onesided(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p)
{
    refine<9>(p,
        [&](const Eigen::Matrix3d &H, auto &JtJ, auto &g, double &cost) {
            const double k = p.distortion_parameter2;
            const double H00 = H(0, 0), H01 = H(0, 1), H02 = H(0, 2);
            const double H10 = H(1, 0), H11 = H(1, 1), H12 = H(1, 2);
            const double H20 = H(2, 0), H21 = H(2, 1);

            for (size_t i = 0; i < x1.size(); i++) {
                const double u = x1[i](0), v = x1[i](1);
                const double p0 = x2[i](0), p1 = x2[i](1);

                const double w = 1.0 / (H20 * u + H21 * v + H(2, 2));
                const double z0 = (H00 * u + H01 * v + H02) * w;
                const double z1 = (H10 * u + H11 * v + H12) * w;

                const double r22 = p0 * p0 + p1 * p1;
                const double df = 1.0 / (1.0 + k * r22);
                const double re = z0 - p0 * df;
                const double ro = z1 - p1 * df;

                const double a = u * w, b = v * w, c = w;
                const double az0 = a * z0, bz0 = b * z0;
                const double az1 = a * z1, bz1 = b * z1;
                accumulate_h_block(JtJ, g, cost,
                    a, b, c, az0, bz0, az1, bz1, re, ro);

                const double r22df2 = r22 * df * df;
                accumulate_dist_col<9, 8>(JtJ, g,
                    a, b, c, az0, bz0, az1, bz1,
                    p0 * r22df2, p1 * r22df2, re, ro);
            }
            return true;
        },
        [&p](const auto &dx) { p.distortion_parameter2 += dx(8); });
}

// Non-linear refinement of transfer error |x2 - pi(H*x1(k))|^2
void refinement_onesided_right(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p)
{
    refine<9>(p,
        [&](const Eigen::Matrix3d &H, auto &JtJ, auto &g, double &cost) {
            const double k = p.distortion_parameter;
            const double H00 = H(0, 0), H01 = H(0, 1), H02 = H(0, 2);
            const double H10 = H(1, 0), H11 = H(1, 1), H12 = H(1, 2);
            const double H20 = H(2, 0), H21 = H(2, 1), H22 = H(2, 2);

            for (size_t i = 0; i < x1.size(); i++) {
                const double u = x1[i](0), v = x1[i](1);
                const double p0 = x2[i](0), p1 = x2[i](1);

                const double r12 = u * u + v * v;
                const double d1 = 1.0 + k * r12;

                const double w = 1.0 / (H20 * u + H21 * v + H22 * d1);
                const double z0 = (H00 * u + H01 * v + H02 * d1) * w;
                const double z1 = (H10 * u + H11 * v + H12 * d1) * w;

                const double re = z0 - p0;
                const double ro = z1 - p1;

                const double a = u * w, b = v * w, c = d1 * w;
                const double az0 = a * z0, bz0 = b * z0;
                const double az1 = a * z1, bz1 = b * z1;
                accumulate_h_block(JtJ, g, cost,
                    a, b, c, az0, bz0, az1, bz1, re, ro);

                const double r12w = r12 * w;
                accumulate_dist_col<9, 8>(JtJ, g,
                    a, b, c, az0, bz0, az1, bz1,
                    r12w * (H02 - z0 * H22), r12w * (H12 - z1 * H22), re, ro);
            }
            return true;
        },
        [&p](const auto &dx) { p.distortion_parameter += dx(8); });
}

void refinement_twosided_equal(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p)
{
    refinement_twosided_impl<true>(x1, x2, p);
}

void refinement_twosided(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p)
{
    refinement_twosided_impl<false>(x1, x2, p);
}

}
