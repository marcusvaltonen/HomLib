#include <Eigen/Dense>
#include <vector>
#include <limits>

#include "PoseLib/misc/sturm.h"

#include "posedata.hpp"
#include "wadenback_common.hpp"


namespace HomLib {
namespace Wadenback3DV2026 {
static inline double unused_equation_residual(
    const double* d,
    double k
);
std::vector<HomLib::PoseData> get_double_sided_equal(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
)
{
    Eigen::Matrix<double, 24, 1> data;
    data << HomLib::Wadenback3DV2026::compute_coeffs(y), HomLib::Wadenback3DV2026::compute_coeffs(x);  // NOTE: y first

    // Compute coefficients
    const double* d = data.data();
    double p[7];
    p[6] = d[1]*d[14] - d[2]*d[13];
    p[5] = d[1]*d[17] - d[2]*d[16] + d[4]*d[14] - d[5]*d[13];
    p[4] = d[1]*d[20] - d[2]*d[19] + d[4]*d[17] - d[5]*d[16] + d[7]*d[14] - d[8]*d[13];
    p[3] = d[1]*d[23] - d[2]*d[22] + d[4]*d[20] - d[5]*d[19] + d[7]*d[17] - d[8]*d[16] + d[10]*d[14] - d[11]*d[13];
    p[2] = d[4]*d[23] - d[5]*d[22] + d[7]*d[20] - d[8]*d[19] + d[10]*d[17] - d[11]*d[16];
    p[1] = d[7]*d[23] - d[8]*d[22] + d[10]*d[20] - d[11]*d[19];
    p[0] = d[10]*d[23] - d[11]*d[22];

    // Find real roots
    double roots[6];
    int nroots = poselib::sturm::bisect_sturm<6>(p, roots);

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
                res = unused_equation_residual(d, k);
                if (res < min_res) {
                    min_res = res;
                    best_id = i;
                }
            }
        }
        
        // Only compute homography for the best one
        double k = roots[best_id];
        std::vector<Eigen::Vector2d> x_undist;
        std::vector<Eigen::Vector2d> y_undist;
        for (int i = 0; i < 4; i++) {
            Eigen::Vector2d tmp;
            tmp = x[i] / (1 + k * x[i].squaredNorm());
            x_undist.push_back(tmp);
            tmp = y[i] / (1 + k * y[i].squaredNorm());
            y_undist.push_back(tmp);
        }
        Eigen::Matrix3d H = HomLib::Wadenback3DV2026::homography_4pt_guo(x_undist, y_undist);

        // Package
        HomLib::PoseData pd;
        pd.homography = H;
        pd.distortion_parameter = k;
        output.push_back(pd);
    } else {
        // Compute homography
        for (int j = 0; j < nroots; j++) {
            HomLib::PoseData pd;
            Eigen::Matrix3d H;
            double k = roots[j];

            std::vector<Eigen::Vector2d> x_undist;
            std::vector<Eigen::Vector2d> y_undist;
            for (int i = 0; i < 4; i++) {
	            Eigen::Vector2d tmp;
	            tmp = x[i] / (1 + k * x[i].squaredNorm());
	            x_undist.push_back(tmp);
	            tmp = y[i] / (1 + k * y[i].squaredNorm());
	            y_undist.push_back(tmp);
            }
            H = HomLib::Wadenback3DV2026::homography_4pt_guo(x_undist, y_undist);

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
    double k
) {
    double k2 = k * k;
    double k3 = k2 * k;
    double k4 = k3 * k;
    double k5 = k4 * k;
    double k6 = k5 * k;

    double res = std::abs(
        (-d[0]*d[14] + d[2]*d[12]) * k6 +
        (-d[0]*d[17] + d[2]*d[15] - d[3]*d[14] + d[5]*d[12]) * k5 +
        (-d[0]*d[20] + d[2]*d[18] - d[3]*d[17] + d[5]*d[15] - d[6]*d[14] + d[8]*d[12]) * k4 +
        (-d[0]*d[23] + d[2]*d[21] - d[3]*d[20] + d[5]*d[18] - d[6]*d[17] + d[8]*d[15] - d[9]*d[14] + d[11]*d[12]) * k3 +
        (-d[3]*d[23] + d[5]*d[21] - d[6]*d[20] + d[8]*d[18] - d[9]*d[17] + d[11]*d[15]) * k2 +
        (-d[6]*d[23] + d[8]*d[21] - d[9]*d[20] + d[11]*d[18]) * k +
        -d[9]*d[23] + d[11]*d[21]
    );
    return res;
}
} // namespace Wadenback3DV2026
} // namespace HomLib
