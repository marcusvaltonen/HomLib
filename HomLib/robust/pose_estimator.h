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

        // xu: xd undistorted with the pose's parameters, yd: yu distorted; false where there is
        // none (see radial.hpp)
        inline bool undistort(const HomLib::PoseData &pose, const Eigen::Vector2d &xd, Eigen::Vector2d &xu) const {
            return HomLib::try_radialundistort(xd, pose.distortion_parameter, xu);
        }
        inline bool distort(const HomLib::PoseData &pose, const Eigen::Vector2d &yu, Eigen::Vector2d &yd) const {
            return HomLib::try_radialdistort(yu, pose.distortion_parameter2, yd);
        }

		// Public, not protected: since C++20, Solver{} for a solver that declares no constructor
		// (an aggregate) initializes this base from outside the solver
		PoseEstimator() = default;
	};
};


template<class Solver>
int HomLib::PoseEstimator<Solver>::estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, std::vector<HomLib::PoseData> *poses) const
{
    return estimate_normalized(x_, y_, poses, [&](const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) {
        static_cast<const Solver*>(this)->solve(x, y, poses, k_range);
    });
}


