#ifndef SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
#define SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
#include <Eigen/Dense>
#include <vector>

namespace HomLib {
namespace Wadenback3DV2026 {
// The first set of (possibly) false roots coming from inv()->adj().
// This computes the determinant of [x{1} x{2} x{3}] (with w_i = 1 + k1*r_i^2), in absolute value.
double compute_residual1(
    const std::vector<Eigen::Vector2d> &x,
    double k1
);
// The second set of (possibly) false roots coming from inv()->adj().
// This computes the determinant of diag(g), i.e. g1*g2*g3, where
// g = adj([x{1} x{2} x{3}]) * x{4} (with w_i = 1 + k1*r_i^2).
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
} // namespace
} // namespace
#endif  // SRC_SOLVERS_WADENBACK_3DV_2026_WADENBACK_COMMON_HPP_
