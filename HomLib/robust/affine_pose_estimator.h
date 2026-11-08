#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "posedata.hpp"
#include "estimator_base.h"

namespace HomLib {

	// We use CRTP here for the solvers.
	template<class Solver>
	class AffinePoseEstimator : public EstimatorBase {
	public:		
		int estimate(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Matrix2d> &A, std::vector<HomLib::PoseData> *poses) const;
		
		inline int minimal_sample_size() const {
			return static_cast<const Solver*>(this)->minimal_sample_size();
		}

		// Public, not protected: since C++20, Solver{} for a solver that declares no constructor
		// (an aggregate) initializes this base from outside the solver
		AffinePoseEstimator() = default;
	};
};


template<class Solver>
int HomLib::AffinePoseEstimator<Solver>::estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, const std::vector<Eigen::Matrix2d> &A_, std::vector<HomLib::PoseData> *poses) const
{
    return estimate_normalized(x_, y_, poses, [&](const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) {
        static_cast<const Solver*>(this)->solve(x, y, A_, poses, k_range);
    });
}


