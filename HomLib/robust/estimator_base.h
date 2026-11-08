#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
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
		// translated: the origin is the distortion center. Returns 0, without changing the points,
		// if f0 is 0 (all of x at the origin) or not finite: the sample has no solution.
		double normalize(std::vector<Eigen::Vector2d> &x, std::vector<Eigen::Vector2d> &y) const {
			if (!normalize_image_coord)
				return 1.0;
			double f0 = 0.0;
			for (size_t i = 0; i < x.size(); i++) {
				f0 += x[i].norm();
			}
			f0 /= x.size();
			f0 /= std::sqrt(2.0);
			if (!(std::isfinite(f0) && f0 > 0.0))
				return 0.0;
			for (size_t i = 0; i < x.size(); i++) {
				x[i] /= f0;
				y[i] /= f0;
			}
			return f0;
		}

		// The estimate() of the three estimators: normalizes copies of x and y, calls
		// solve(x, y) on them (which adds the solutions to poses), and denormalizes the solutions;
		// their number. No solution if the sample can't be normalized.
		template <class Solve>
		int estimate_normalized(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_,
		                        std::vector<HomLib::PoseData> *poses, const Solve &solve) const {
			std::vector<Eigen::Vector2d> x = x_;
			std::vector<Eigen::Vector2d> y = y_;
			const double f0 = normalize(x, y);
			poses->clear();
			if (f0 == 0.0)
				return 0;
			solve(x, y);
			denormalize(poses, f0);
			return static_cast<int>(poses->size());
		}

		// The solutions back in the input's coordinates. The solvers return only finite solutions
		// (tested for every solver in test_finite_solutions.cpp); checked here in Debug builds.
		void denormalize(std::vector<HomLib::PoseData> *poses, double f0) const {
			assert(std::all_of(poses->begin(), poses->end(), [](const HomLib::PoseData &p) {
				return p.homography.allFinite() && std::isfinite(p.distortion_parameter) &&
				       std::isfinite(p.distortion_parameter2) && std::isfinite(p.focal_length);
			}));
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
		}
	};
};
