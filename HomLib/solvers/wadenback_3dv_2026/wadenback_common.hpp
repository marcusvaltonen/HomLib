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
// Finds the real roots of p[0] + p[1] x + ... + p[N] x^N, writes them to roots (at most N), and
// returns their number. It calls PoseLib's bisect_sturm on q(z) = p(s (z + c)).
//
// bisect_sturm stops on absolute thresholds. Ridders' method stops when the bracket is below 1e-3
// or after 30 iterations. Newton's method stops when |p(x)| / |p_N| < tol (1e-10 by default). The
// solvers' polynomials have small roots, close roots and large spurious roots. bisect_sturm finds
// the small and close roots coarsely, and the large ones make its search interval huge. Three
// changes fix this:
// - s = max(|p[0] / p[N]|^(1/N), |p[1] / p[N]|^(1/(N-1))) brings the typical root near 1.
// - tol = 1e-14 lets bisect_sturm's Newton steps reach the rounding level of q. This also makes
//   the close roots accurate.
// - Without distortion, p has a root at x = 0, which is bisect_sturm's first bisection point.
//   There, p = 0 up to rounding, and bisect_sturm can lose the root. The offset c moves the root
//   to z = -c.
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
    // Any small offset works, and 0.0123 has no special meaning. The offset only has to move a root
    // at x = 0 off bisect_sturm's bisection points z = B n / 2^j. z = -c is one of them only if
    // c 2^j / B is an integer.
    const double c = 0.0123;
    for (int i = 0; i < N; i++)   // Taylor shift: q(z) becomes q(z + c)
        for (int j = N - 1; j >= i; j--)
            q[j] += c * q[j + 1];
    int nroots = poselib::sturm::bisect_sturm<N>(q, roots, 1e-14);
    for (int i = 0; i < nroots; i++)
        roots[i] = (roots[i] + c) * s;
    return nroots;
}
// Guo's method (2019) for computing a homography from four points using an explicit method.
// False, H unchanged, if three of the first four points of x are collinear (H can still be
// non-finite, on overflow).
bool homography_4pt_guo(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    Eigen::Matrix3d &H
);
// The same, from the first four points undistorted by the division model, with k1 for x and k2
// for y; false also if a point can't be undistorted (1 + k |x|^2 = 0), or H isn't finite.
bool homography_4pt_guo_undistorted(
    const std::vector<Eigen::Vector2d> &x, double k1,
    const std::vector<Eigen::Vector2d> &y, double k2,
    Eigen::Matrix3d &H
);
} // namespace
} // namespace
#endif  // SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
