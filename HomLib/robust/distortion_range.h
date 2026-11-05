#pragma once

#include <limits>
#include <vector>
#include "posedata.hpp"

namespace HomLib {

// Range check on the distortion parameters of a solver's solutions, in the coordinates of the
// points it gets: both parameters must be in [k_min, k_max]. A parameter equal to 0, the one that
// a one-sided solver doesn't estimate, passes. DistortionRange() checks nothing (the solvers'
// default argument); the robust estimators use [-1, 0.1] by default.
//
// The robust estimators hold the range, k_range (EstimatorBase), and normalize each sample before
// calling the solver: they divide both images' points by f0 = mean(|x|)/sqrt(2), computed from the
// first image's points x (measured from the distortion center, the origin). The range is in these
// normalized coordinates, whatever the units of the input. The solvers check it before computing
// the homographies of the rejected roots, so RANSAC has fewer hypotheses to score. It only
// selects the hypotheses: the refinement isn't constrained.
//
// The default [-1, 0.1] is safe: to first order, the relative distortion at the image corners
// is D = k r_c^2, and for points spread over an image, r_c is about 2.6 in the normalized
// coordinates, so the range only rejects more than about 70% of pincushion distortion there,
// whatever the focal length (see the solver guide).
//
// TODO: Kukelova et al.'s 5-point solvers normalize their points again (with the larger of the
// two images' scales; they check the range after scaling back): normalize once, in one place.
struct DistortionRange {
    bool enabled = false;
    double k_min = -std::numeric_limits<double>::infinity();
    double k_max = std::numeric_limits<double>::infinity();

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
