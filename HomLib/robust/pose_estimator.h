#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "posedata.hpp"
#include "estimator_base.h"
#include "solver_id.h"
#include "radial.hpp"
#include "refinement.hpp"

namespace HomLib {

	// We use CRTP here for the solvers.
	template<class Solver>
	class PoseEstimator : public EstimatorBase {
	public:
		static constexpr bool supports_overdetermined = false;
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

	template <class Solver>
	class NoDistortionPointEstimator : public PoseEstimator<Solver> {
	public:
		static constexpr DistortionCase DISTORTION_CASE = NO_DISTORTION;

		int minimal_sample_size() const {
			return Solver::minimal_sample_size_value;
		}

		inline bool undistort(const PoseData &, const Eigen::Vector2d &xd, Eigen::Vector2d &xu) const {
			xu = xd;
			return true;
		}

		inline bool distort(const PoseData &, const Eigen::Vector2d &yu, Eigen::Vector2d &yd) const {
			yd = yu;
			return true;
		}

		inline void refine(PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
			refinement_no_dist(x, y, pose);
		}
	};

	template <class Solver>
	class OneSidedLeftPointEstimator : public PoseEstimator<Solver> {
	public:
		static constexpr DistortionCase DISTORTION_CASE = ONE_SIDED_LEFT;

		int minimal_sample_size() const {
			return Solver::minimal_sample_size_value;
		}

		inline bool undistort(const PoseData &, const Eigen::Vector2d &xd, Eigen::Vector2d &xu) const {
			xu = xd;
			return true;
		}

		inline void refine(PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
			refinement_onesided(x, y, pose);
		}
	};

	template <class Solver>
	class TwoSidedEqualPointEstimator : public PoseEstimator<Solver> {
	public:
		static constexpr DistortionCase DISTORTION_CASE = TWO_SIDED_EQUAL;

		int minimal_sample_size() const {
			return Solver::minimal_sample_size_value;
		}

		inline void refine(PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
			refinement_twosided_equal(x, y, pose);
		}
	};

	template <class Solver>
	class TwoSidedPointEstimator : public PoseEstimator<Solver> {
	public:
		static constexpr DistortionCase DISTORTION_CASE = TWO_SIDED;

		int minimal_sample_size() const {
			return Solver::minimal_sample_size_value;
		}

		inline void refine(PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
			refinement_twosided(x, y, pose);
		}
	};
};


template<class Solver>
int HomLib::PoseEstimator<Solver>::estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, std::vector<HomLib::PoseData> *poses) const
{
    return estimate_normalized(x_, y_, poses, [&](const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) {
        static_cast<const Solver*>(this)->solve(x, y, poses, k_range);
    });
}


