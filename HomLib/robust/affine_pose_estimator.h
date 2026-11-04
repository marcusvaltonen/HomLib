#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "posedata.hpp"
#include "distortion_range.h"

namespace HomLib {

	// We use CRTP here for the solvers.
	template<class Solver>
	class AffinePoseEstimator {
	public:		
		int estimate(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses) const;
		
		inline int minimal_sample_size() const {
			return static_cast<const Solver*>(this)->minimal_sample_size();
		}

		// Options
		bool normalize_image_coord = true;

		// Optional range check on the distortion parameters of the solutions (off by default):
		// keeps those with both parameters in [k_min, k_max] (a parameter equal to 0, not estimated
		// by one-sided solvers, passes). It isn't part of the solvers' methods, but it removes
		// implausible solutions, so RANSAC has fewer hypotheses to score. The range applies in the
		// normalized coordinates (each sample is divided by mean(|x|)/sqrt(2), with
		// normalize_image_coord), or with k_range_in_input_units, in the units of the points given to
		// estimate(). It only selects the hypotheses: refinement isn't constrained. The range is
		// passed to the solver's solve(), which may apply it to each root before computing its
		// homography (and, if it keeps one root of several, before choosing it); the solutions of
		// all the solvers are then filtered.
		bool use_k_range = false;
		double k_min = -1.0;
		double k_max = 0.1;
		bool k_range_in_input_units = false;

	protected:
		AffinePoseEstimator() = default;
	};
};


template<class Solver>
int HomLib::AffinePoseEstimator<Solver>::estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, const std::vector<Eigen::Matrix2d> &A_, std::vector<HomLib::PoseData> *poses) const
{
    std::vector<Eigen::Vector2d> x = x_;
    std::vector<Eigen::Vector2d> y = y_;
    std::vector<Eigen::Matrix2d> A = A_;

    // Rescale image plane
    double f0 = 0.0;
    if (normalize_image_coord) {
        // TODO: Consider full Hartley normalization, i.e. also translate.
        for (size_t i = 0; i < x.size(); i++) {
            f0 += x[i].norm();
        }
        f0 /= x.size();
        f0 /= std::sqrt(2.0);
        for (size_t i = 0; i < x.size(); i++) {
            x[i] /= f0;
            y[i] /= f0;
        }
    }

    // Call solver implementation
    poses->clear();
    // The range in the solver's (normalized) coordinates
    HomLib::DistortionRange range;
    if (use_k_range) {
        const double s = (k_range_in_input_units && normalize_image_coord) ? f0 * f0 : 1.0;
        range = HomLib::DistortionRange(k_min * s, k_max * s);
    }
    int n_sols = static_cast<const Solver*>(this)->solve(x, y, A, poses, range);
    if (use_k_range && !k_range_in_input_units) {
        HomLib::filter_distortion_range(poses, HomLib::DistortionRange(k_min, k_max));
        n_sols = poses->size();
    }

    // Revert image coordinate scaling
    if (normalize_image_coord) {
        double f02 = f0 * f0;
        for (size_t i = 0; i < poses->size(); ++i) {
            (*poses)[i].homography(0,0) *= f0;
            (*poses)[i].homography(0,1) *= f0;
            (*poses)[i].homography(0,2) *= f02;
            (*poses)[i].homography(1,0) *= f0;
            (*poses)[i].homography(1,1) *= f0;
            (*poses)[i].homography(1,2) *= f02;
            (*poses)[i].homography(2,2) *= f0;
            (*poses)[i].distortion_parameter /= f02;
            (*poses)[i].distortion_parameter2 /= f02;
        }
    }

    if (use_k_range && k_range_in_input_units) {
        HomLib::filter_distortion_range(poses, HomLib::DistortionRange(k_min, k_max));
        n_sols = poses->size();
    }

    return n_sols;
}


