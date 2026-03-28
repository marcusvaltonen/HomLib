#include <Eigen/Dense>
#include <vector>
#include <limits>

#include "PoseLib/misc/univariate.h"

#include "posedata.hpp"
#include "wadenback_common.hpp"


namespace HomLib {
namespace Wadenback2025 {
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
    Eigen::Matrix<double, 12, 1> data = HomLib::Wadenback2025::compute_coeffs(y);
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

    // Only check residual if nroots > 1
    double min_res = std::numeric_limits<double>::max();
    double res;
    int best_id = 0;
    std::vector<HomLib::PoseData> output;


    if (nroots == 0) {
        return output;
    }
    
    if (extra_check) {
        if (nroots > 1) {
            for (int i = 0; i < nroots; i++) {
                double k = roots[i];
                res = unused_equation_residual(d, rhs, k);
                if (res < min_res) {
                    min_res = res;
                    best_id = i;
                }
            }
        }
        double k = roots[best_id];
        std::vector<Eigen::Vector2d> y_undist;
        for (int i = 0; i < 4; i++) {
            Eigen::Vector2d tmp = y[i] / (1 + k * y[i].squaredNorm());
            y_undist.push_back(tmp);
        }
        Eigen::Matrix3d H = HomLib::Wadenback2025::homography_4pt_guo(x, y_undist);
        // Package
        HomLib::PoseData pd;
        pd.homography = H;
        pd.distortion_parameter = k;

        output.push_back(pd);
    } else {
        for (int j = 0; j < nroots; j++) {
            // Compute homography
            HomLib::PoseData pd;
            Eigen::Matrix3d H;
            double k = roots[j];

            std::vector<Eigen::Vector2d> y_undist;
            for (int i = 0; i < 4; i++) {
	            Eigen::Vector2d tmp = y[i] / (1 + k * y[i].squaredNorm());
	            y_undist.push_back(tmp);
            }
            H = HomLib::Wadenback2025::homography_4pt_guo(x, y_undist);
            // Package
            pd.homography = H;
            pd.distortion_parameter = k;

            output.push_back(pd);
        }
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

} // namespace Wadenback2025
} // namespace HomLib
