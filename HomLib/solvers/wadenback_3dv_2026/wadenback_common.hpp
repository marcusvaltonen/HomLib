#ifndef SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
#define SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <vector>
#include "PoseLib/misc/sturm.h"

namespace HomLib {
namespace Wadenback3DV2026 {
// The first set of (possibly) false roots coming from inv()->adj(), where det([x{1} x{2} x{3}])
// vanishes (with w_i = 1 + k1*r_i^2). This computes |det([x{1} x{2} x{3}])| divided by the product
// of the row norms: a scale-free value in [0, 1], 0 at these roots.
double compute_residual1(
    const std::vector<Eigen::Vector2d> &x,
    double k1
);
// The second set of (possibly) false roots coming from inv()->adj(), where the determinant of
// diag(g) vanishes, i.e. a component of g = adj([x{1} x{2} x{3}]) * x{4} (with w_i = 1 + k1*r_i^2).
// This computes min_i |g_i| / |g|: a scale-free value in [0, 1], 0 at these roots.
double compute_residual2(
    const std::vector<Eigen::Vector2d> &x,
    double k1
);
// Coefficients (row-major order) of the 3x4 matrix N writing LHS/RHS as N * [k^3; k^2; k; 1]
Eigen::Matrix<double, 12, 1> compute_coeffs(
    const std::vector<Eigen::Vector2d> &x
);
// Writes the real roots of p[0] + p[1] x + ... + p[N] x^N to roots (at most N) and returns their
// number. Uses PoseLib's bisect_sturm on q(z) = p(s (z + c)), then one Newton step on p per root.
//
// bisect_sturm searches [-B, B], B = 1 + max |p_i / p_N|, and its stops are absolute (Ridders'
// method until the bracket is below 1e-3 or after 30 iterations, then Newton until
// |p(x)| / |p_N| < 1e-10). The solvers' polynomials have small, plausible roots and large,
// spurious ones: these make B huge, so the iterations can stop far from a root; small roots are
// found coarsely; and a root at 0 (no distortion) is the first bisection point, where the Sturm
// count (p from the Sturm sequence) and Ridders' method (p from polyval) see p = 0 up to rounding,
// possibly with different signs: then the root is counted in a half where Ridders' method finds
// no sign change, and is dropped. s = max(|p[0] / p[N]|^(1/N), |p[1] / p[N]|^(1/(N-1))) brings the
// typical root magnitude near 1, and the small offset c = 0.0123 (any small value works) moves a
// root at x = 0 to z = -c, a bisection point (B n / 2^j) only if c 2^j / B is an integer.
// The Newton step, kept if it reduces |p|, restores the accuracy lost when s > 1 (on q, the stop
// on |p| is looser by s^N). It is also kept only if it moves x by at most 1e-3 s, the size of
// bisect_sturm's last bracket: near a multiple root or two close roots, p' is close to 0 and the
// step can jump to another root, which would then be returned twice and this one lost.
template <int N>
int bisect_sturm_scaled(const double *p, double *roots) {
    static_assert(N >= 2, "degree at least 2");
    double s = std::max(std::pow(std::abs(p[0] / p[N]), 1.0 / N),
                        std::pow(std::abs(p[1] / p[N]), 1.0 / (N - 1)));
    if (!(s > 0) || !std::isfinite(s))
        s = 1;
    double q[N + 1];
    double sk = 1;
    for (int i = 0; i <= N; i++) {
        q[i] = p[i] * sk;
        sk *= s;
    }
    const double c = 0.0123;
    for (int i = 0; i < N; i++)   // Taylor shift: q(z) becomes q(z + c)
        for (int j = N - 1; j >= i; j--)
            q[j] += c * q[j + 1];
    int nroots = poselib::sturm::bisect_sturm<N>(q, roots);
    for (int i = 0; i < nroots; i++) {
        const double x = (roots[i] + c) * s;
        double f = p[N], df = 0;   // Horner: f = p(x), df = p'(x)
        for (int j = N - 1; j >= 0; j--) {
            df = df * x + f;
            f = f * x + p[j];
        }
        roots[i] = x;
        if (df != 0) {
            const double x1 = x - f / df;
            double f1 = p[N];
            for (int j = N - 1; j >= 0; j--)
                f1 = f1 * x1 + p[j];
            if (std::abs(f1) < std::abs(f) && std::abs(x1 - x) <= 1e-3 * s)
                roots[i] = x1;
        }
    }
    return nroots;
}
// Guo's method (2019) for computing a homography from four points using an explicit method.
Eigen::Matrix3d homography_4pt_guo(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y
);
} // namespace
} // namespace
#endif  // SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
