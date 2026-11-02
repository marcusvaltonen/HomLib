#pragma once
#ifdef HOMLIB_HAS_SUPERANSAC

#include <Eigen/Dense>
#include "posedata.hpp"
#include "models/model.h"
#include "utils/types.h"

namespace HomLib {

inline superansac::models::Model pose_to_model(const PoseData& pose) {
    superansac::models::Model model;
    auto& data = model.getMutableData();
    data.resize(4, 4);
    data.setZero();
    data.block<3, 3>(0, 0) = pose.homography;
    data(3, 0) = pose.focal_length;
    data(3, 1) = pose.distortion_parameter;
    data(3, 2) = pose.distortion_parameter2;
    return model;
}

inline PoseData model_to_pose(const superansac::models::Model& model) {
    const auto& data = model.getData();
    PoseData pose;
    pose.homography = data.block<3, 3>(0, 0);
    pose.focal_length = data(3, 0);
    pose.distortion_parameter = data(3, 1);
    pose.distortion_parameter2 = data(3, 2);
    return pose;
}

} // namespace HomLib

#endif // HOMLIB_HAS_SUPERANSAC
