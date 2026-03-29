#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "posedata.hpp"

namespace HomLib {

	// We use CRTP here for the solvers.
	template<class Solver>
	class PoseEstimator {
	public:		
		int estimate(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const;
		
		inline int minimal_sample_size() const {
			return static_cast<const Solver*>(this)->minimal_sample_size();
		}

		// Options
		bool normalize_image_coord = true;

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
    int n_sols = static_cast<const Solver*>(this)->solve(x, y, poses);

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

    return n_sols;
}


