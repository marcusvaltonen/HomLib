#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include <limits>

#include "PoseLib/misc/univariate.h"

#include "posedata.hpp"
#include "wadenback_common.hpp"


namespace HomLib {
namespace Wadenback3DV2026 {
static inline double unused_equation_residual(
    const double* d,
    const Eigen::Vector3d &rhs,
    double k
);
static inline Eigen::Vector3d get_rhs(
    const std::vector<Eigen::Vector2d> &x
);
std::vector<HomLib::PoseData> get_one_sided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
)
{
    Eigen::Matrix<double, 12, 1> data = HomLib::Wadenback3DV2026::compute_coeffs(y);
    Eigen::Vector3d rhs = get_rhs(x);

    // Compute coefficients
    // TODO: Plenty of coefficients that are not necessary to compute here... (or maybe used in "unused eq later..")
    const double* d = data.data();
    const double c3 = -rhs(1)*d[2] + rhs(2)*d[1];
    const double c2 = -rhs(1)*d[5] + rhs(2)*d[4];
    const double c1 = -rhs(1)*d[8] + rhs(2)*d[7];
    const double c0 = -rhs(1)*d[11] + rhs(2)*d[10];

    // Find real roots
    double inv_c3 = 1.0 / c3;
    double roots[3];
    int nroots = poselib::univariate::solve_cubic_real(c2 * inv_c3, c1 * inv_c3, c0 * inv_c3, roots);

    std::vector<HomLib::PoseData> output;


    if (nroots == 0) {
        return output;
    }
    
    // With extra_check, only the root with the smallest residual among those that a homography
    // can go with (see homography_of_root): try the roots by increasing residual
    int order[3];
    for (int i = 0; i < nroots; i++)
        order[i] = i;
    if (extra_check && nroots > 1) {
        double res_of[3];
        for (int i = 0; i < nroots; i++)
            res_of[i] = unused_equation_residual(d, rhs, roots[i]);
        std::sort(order, order + nroots, [&res_of](int a, int b) { return res_of[a] < res_of[b]; });
    }
    for (int j = 0; j < nroots; j++) {
        const double k = roots[order[j]];
        HomLib::PoseData pd;
        if (!HomLib::Wadenback3DV2026::homography_of_root(x, y, 0.0, k, &pd.homography))
            continue;   // a false root
        pd.distortion_parameter2 = k;
        output.push_back(pd);
        if (extra_check)
            break;
    }

	return output;
}

static inline double unused_equation_residual(
    const double* d,
    const Eigen::Vector3d &rhs,
    double k
) {
    double k2 = k * k;
    double k3 = k2 * k;
    double res = std::abs(
        (rhs(0)*d[2] - rhs(2)*d[0]) * k3 +
        (rhs(0)*d[5] - rhs(2)*d[3]) * k2 +
        (rhs(0)*d[8] - rhs(2)*d[6]) * k +
        rhs(0)*d[11] - rhs(2)*d[9]
    );
    return res;
}

static inline Eigen::Vector3d get_rhs(
    const std::vector<Eigen::Vector2d> &x
) {
    double m11 = x[1][1] - x[2][1];
    double m12 = -x[1][0] + x[2][0];
    double m13 = x[1][0]*x[2][1] - x[1][1]*x[2][0];
    double m21 = -x[0][1] + x[2][1];
    double m22 = x[0][0] - x[2][0];
    double m23 = -x[0][0]*x[2][1] + x[0][1]*x[2][0];
    double m31 = x[0][1] - x[1][1];
    double m32 = -x[0][0] + x[1][0];
    double m33 = x[0][0]*x[1][1] - x[0][1]*x[1][0];
    double g1 = m13 + m11*x[3][0] + m12*x[3][1];
    double g2 = m23 + m21*x[3][0] + m22*x[3][1];
    double g3 = m33 + m31*x[3][0] + m32*x[3][1];
    double h1 = m13 + m11*x[4][0] + m12*x[4][1];
    double h2 = m23 + m21*x[4][0] + m22*x[4][1];
    double h3 = m33 + m31*x[4][0] + m32*x[4][1];

    Eigen::Vector3d rhs;
    rhs << g2*g3*h1,
           g1*g3*h2,
           g1*g2*h3;

    return rhs;
}

} // namespace Wadenback3DV2026
} // namespace HomLib
