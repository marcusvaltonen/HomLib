#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "posedata.hpp"
#include "estimator_base.h"

namespace HomLib {

	// We use CRTP here for the solvers.
	template<class Solver>
	class PoseEstimator : public EstimatorBase {
	public:		
		int estimate(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const;
		
		inline int minimal_sample_size() const {
			return static_cast<const Solver*>(this)->minimal_sample_size();
		}

		// Distortion functions using the pose's parameters
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter2);
            return yd;
        }

	protected:
		PoseEstimator() = default;
	};
};


template<class Solver>
int HomLib::PoseEstimator<Solver>::estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, std::vector<HomLib::PoseData> *poses) const
{
    std::vector<Eigen::Vector2d> x = x_;
    std::vector<Eigen::Vector2d> y = y_;
    const double f0 = normalize(x, y);
    poses->clear();
    int n_sols = static_cast<const Solver*>(this)->solve(x, y, poses, k_range);
    denormalize(poses, f0);
    return n_sols;
}


