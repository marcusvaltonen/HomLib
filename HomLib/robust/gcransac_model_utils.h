#pragma once
#ifdef HOMLIB_HAS_GCRANSAC

#include <vector>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>
#include "posedata.hpp"
#include "model.h"

namespace HomLib {

// ---------------------------------------------------------------------------
// PoseData <-> gcransac::Model conversion
// ---------------------------------------------------------------------------

inline gcransac::Model pose_to_gc_model(const PoseData& pose) {
    gcransac::Model model;
    model.descriptor.resize(4, 4);
    model.descriptor.setZero();
    model.descriptor.block<3, 3>(0, 0) = pose.homography;
    model.descriptor(3, 0) = pose.focal_length;
    model.descriptor(3, 1) = pose.distortion_parameter;
    model.descriptor(3, 2) = pose.distortion_parameter2;
    return model;
}

inline PoseData gc_model_to_pose(const gcransac::Model& model) {
    PoseData pose;
    pose.homography = model.descriptor.block<3, 3>(0, 0);
    pose.focal_length = model.descriptor(3, 0);
    pose.distortion_parameter = model.descriptor(3, 1);
    pose.distortion_parameter2 = model.descriptor(3, 2);
    return pose;
}

} // namespace HomLib

#endif // HOMLIB_HAS_GCRANSAC
