#pragma once

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"

namespace HomLib {

template<class Solver>
class RansacEstimator {
public:
	RansacEstimator(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, Solver est) {
		x = x_;
		y = y_;
		solver = est;
	}

	inline int min_sample_size() const {
		return solver.minimal_sample_size();
	}
	inline int non_minimal_sample_size() const {
		return solver.minimal_sample_size() * 2;
	}
	inline int num_data() const {
		return x.size();
	}

	int MinimalSolver(const std::vector<int>& sample,
		std::vector<HomLib::PoseData>* poses) const {
		
		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;

		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
		}
		solver.estimate(xx, yy, poses);
				
		return poses->size();
	}

	// Returns 0 if no model could be estimated and 1 otherwise.
	int NonMinimalSolver(const std::vector<int>& sample, HomLib::PoseData* pose) const {
		if (!use_non_minimal)
			return 0;

		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;	

		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
		}

		// Call minimal solver
		std::vector<HomLib::PoseData> poses;
		std::vector<Eigen::Vector2d> xs;
		std::vector<Eigen::Vector2d> ys;
		for (int i = 0; i < min_sample_size(); i++) {
			xs.push_back(xx[i]);
			ys.push_back(yy[i]);
		}
		solver.estimate(xs, ys, &poses);

		// for all pose candidates compute score
		double best_score = std::numeric_limits<double>::max();
		int best_idx = -1;

		for (size_t i = 0; i < poses.size(); ++i) {
			double score = 0;
			for (size_t j = 0; j < sample.size(); ++j)
				score += EvaluateModelOnPoint(poses[i], sample[j]);
			if (score < best_score) {
				best_score = score;
				best_idx = i;
			}
		}

		if (best_idx != -1) {
			*pose = poses[best_idx];
			return 1;
		} else {
			return 0;
		}
	}

	// Evaluates the line on the i-th data point.
	double EvaluateModelOnPoint(const HomLib::PoseData& pose, int i) const {
		// Rectify
		Eigen::Vector2d z = solver.undistort(pose, x[i]);
	
		// Compute reprojection error
		Eigen::Vector3d Z = pose.homography * z.homogeneous();
		z = Z.hnormalized();
		z = solver.distort(pose, z);
		
		// std::cout << "evaluate " << i << ": " << (y[i] - z).squaredNorm() << std::endl;

		return (y[i] - z).squaredNorm();
	}

	// Linear least squares solver. Calls NonMinimalSolver.
	inline void LeastSquares(const std::vector<int>& sample, HomLib::PoseData* p) const {
		if (!use_local_opt)
			return;
		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;	

		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
		}
		solver.refine(*p, xx, yy);
	}


	bool use_non_minimal = true;
	bool use_local_opt = true;
private:
	Solver solver;
	std::vector<Eigen::Vector2d> x;
	std::vector<Eigen::Vector2d> y;

};

}
