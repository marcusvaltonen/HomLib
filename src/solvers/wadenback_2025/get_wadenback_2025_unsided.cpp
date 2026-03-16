
#include <Eigen/Dense>
#include "posedata.hpp"
#include "wadenback_common.hpp"
#include "pose_estimator.h"

namespace HomLib {
namespace Wadenback2025 {
std::vector<HomLib::PoseData> get_unsided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
)
{
    Eigen::Matrix3d H = HomLib::Wadenback2025::homography_4pt_guo(x, y);

    HomLib::PoseData pd;
    pd.homography = H;
    return {pd};
}
}
}
