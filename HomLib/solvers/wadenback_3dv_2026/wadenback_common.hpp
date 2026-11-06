#ifndef SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
#define SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
#include <Eigen/Dense>
#include <vector>

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
// Guo's method (2019) for computing a homography from four points using an explicit method.
Eigen::Matrix3d homography_4pt_guo(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y
);
// The homography of a root, from the first 5 points of x and y undistorted with k1 and k2: Guo's
// method on the 4 of them in the best general position (see general_position.hpp). Returns false
// for a root that can't be a solution: when a point has 1 + k r^2 <= 0 (no undistorted point),
// when the 4 points don't determine a homography with every point in front of both cameras, or
// when H isn't finite.
bool homography_of_root(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    double k1,
    double k2,
    Eigen::Matrix3d *H
);
} // namespace
} // namespace
#endif  // SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
