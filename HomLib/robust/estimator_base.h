#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include "posedata.hpp"
#include "distortion_range.h"
#include "solver_id.h"

namespace HomLib {

	// What the robust estimators (PoseEstimator, AffinePoseEstimator, OrientationPoseEstimator)
	// share: the normalization of each sample, and the range check on the distortion, in the
	// normalized coordinates (see DistortionRange).
	class EstimatorBase {
	public:
		// Options
		bool normalize_image_coord = true;
		HomLib::DistortionRange k_range = HomLib::DistortionRange(-1.0, 0.1);   // in the input's coordinates if normalize_image_coord is false

	protected:
		// With normalize_image_coord, divides x and y by f0 = mean(|x|)/sqrt(2), so that the mean
		// distance of x to the origin is sqrt(2), and returns f0 (otherwise 1). The points aren't
		// translated: the origin is the distortion center.
		double normalize(std::vector<Eigen::Vector2d> &x, std::vector<Eigen::Vector2d> &y) const {
			if (!normalize_image_coord)
				return 1.0;
			double f0 = 0.0;
			for (size_t i = 0; i < x.size(); i++) {
				f0 += x[i].norm();
			}
			f0 /= x.size();
			f0 /= std::sqrt(2.0);
			for (size_t i = 0; i < x.size(); i++) {
				x[i] /= f0;
				y[i] /= f0;
			}
			return f0;
		}

		// The solutions back in the input's coordinates
		void denormalize(std::vector<HomLib::PoseData> *poses, double f0) const {
			if (!normalize_image_coord)
				return;
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
	};
};
