#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "posedata.hpp"

namespace HomLib {

	// We use CRTP here for the solvers.
	template<class Solver>
	class OrientationPoseEstimator {
	public:		
		int estimate(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, const std::vector<Eigen::Vector2d> &ori, std::vector<HomLib::PoseData> *poses) const;
		
		inline int minimal_sample_size() const {
			return static_cast<const Solver*>(this)->minimal_sample_size();
		}

		// Options
		bool normalize_image_coord = true;

	protected:
		OrientationPoseEstimator() = default;
	};
};


template<class Solver>
int HomLib::OrientationPoseEstimator<Solver>::estimate(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, const std::vector<Eigen::Vector2d> &ori_, std::vector<HomLib::PoseData> *poses) const
{
    std::vector<Eigen::Vector2d> x = x_;
    std::vector<Eigen::Vector2d> y = y_;
    std::vector<Eigen::Vector2d> ori = ori_;

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
    int n_sols = static_cast<const Solver*>(this)->solve(x, y, ori, poses);

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


