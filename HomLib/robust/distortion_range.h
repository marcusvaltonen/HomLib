#pragma once

#include <vector>
#include "posedata.hpp"

namespace HomLib {

// Optional range check on the distortion parameters, the same for the solvers and the
// estimators (PoseEstimator, AffinePoseEstimator, OrientationPoseEstimator): the parameters
// must both be in [k_min, k_max]. A parameter equal to 0, the one that a one-sided solver
// doesn't estimate, passes. Disabled by default (then k_min and k_max aren't used).
struct DistortionRange {
    bool enabled = false;
    double k_min = -1.0;
    double k_max = 0.1;

    DistortionRange() = default;
    DistortionRange(double k_min_, double k_max_) : enabled(true), k_min(k_min_), k_max(k_max_) {}

    bool contains(double k) const
    {
        return !enabled || k == 0 || (k >= k_min && k <= k_max);
    }
    bool contains(double k1, double k2) const
    {
        return contains(k1) && contains(k2);
    }
};

// Keeps the solutions whose distortion parameters are in the range
inline void filter_distortion_range(std::vector<HomLib::PoseData> *poses, const DistortionRange &range)
{
    std::vector<HomLib::PoseData> kept;
    for (size_t i = 0; i < poses->size(); i++)
        if (range.contains((*poses)[i].distortion_parameter, (*poses)[i].distortion_parameter2))
            kept.push_back((*poses)[i]);
    poses->swap(kept);
}

}  // namespace HomLib
