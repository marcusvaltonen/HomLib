#include <Eigen/Dense>
#include <algorithm>
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

    std::vector<HomLib::PoseData> output;

    if (nroots == 0) {
        return output;
    }
    
    // With extra_check, only the root with the smallest residual among those that a homography
    // can go with (see homography_of_root): try the roots by increasing residual
    int order[6];
    for (int i = 0; i < nroots; i++)
        order[i] = i;
    if (extra_check && nroots > 1) {
        double res_of[6];
        for (int i = 0; i < nroots; i++)
            res_of[i] = unused_equation_residual(d, roots[i]);
        std::sort(order, order + nroots, [&res_of](int a, int b) { return res_of[a] < res_of[b]; });
    }
    for (int j = 0; j < nroots; j++) {
        const double k = roots[order[j]];
        HomLib::PoseData pd;
        if (!HomLib::Wadenback3DV2026::homography_of_root(x, y, k, k, &pd.homography))
            continue;   // a false root
        pd.distortion_parameter = k;
        pd.distortion_parameter2 = k;
        output.push_back(pd);
        if (extra_check)
            break;
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
