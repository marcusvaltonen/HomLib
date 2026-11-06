#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "posedata.hpp"
#include "estimator_base.h"

namespace HomLib {

	// We use CRTP here for the solvers.
	template<class Solver>
	class OrientationPoseEstimator : public EstimatorBase {
	public:		
		int estimate(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Vector2d> &ori, std::vector<HomLib::PoseData> *poses) const;
		
		inline int minimal_sample_size() const {
			return static_cast<const Solver*>(this)->minimal_sample_size();
		}

	protected:
		OrientationPoseEstimator() = default;
	};
};


template<class Solver>
int HomLib::OrientationPoseEstimator<Solver>::estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, const std::vector<Eigen::Vector2d> &ori_, std::vector<HomLib::PoseData> *poses) const
{
    std::vector<Eigen::Vector2d> x = x_;
    std::vector<Eigen::Vector2d> y = y_;
    const double f0 = normalize(x, y);
    poses->clear();
    int n_sols = static_cast<const Solver*>(this)->solve(x, y, ori_, poses, k_range);
    denormalize(poses, f0);
    return n_sols;
}


