#pragma once

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "solver_utils.h"
namespace HomLib {

template<class Solver>
class OrientationRansacEstimator {
public:
	OrientationRansacEstimator(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, const std::vector<Eigen::Vector2d> &ori_, Solver est) {
		x = x_;
		y = y_;
		ori = ori_;
		solver = est;
	}

	inline int min_sample_size() const {
		return solver.minimal_sample_size();
	}
	inline int non_minimal_sample_size() const {
		return solver.minimal_sample_size() * 4;
	}
	inline int num_data() const {
		return x.size();
	}

	int MinimalSolver(const std::vector<int>& sample,
		std::vector<HomLib::PoseData>* poses) const {
		
		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;
		std::vector<Eigen::Vector2d> oo;

		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
			oo.push_back(ori[sample[i]]);
		}
		solver.estimate(xx, yy, oo, poses);
		
		//std::cout << "H =\n" << (*poses)[0].homography << std::endl;
				
		return poses->size();
	}

	int NonMinimalSolver(const std::vector<int>& sample, HomLib::PoseData* pose) const {
		if (!use_non_minimal)
			return 0;

		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;
		std::vector<Eigen::Vector2d> oo;
		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
			oo.push_back(ori[sample[i]]);
		}

		auto est_fn = [&](std::vector<HomLib::PoseData>* poses) {
			solver.estimate(xx, yy, oo, poses);
		};
		if (!detail::lomsac_non_minimal_full(solver, xx, yy, est_fn, *pose))
			return 0;
		solver.refine(*pose, xx, yy);
		return 1;
	}

	double EvaluateModelOnPoint(const HomLib::PoseData& pose, int i) const {
		return detail::compute_squared_residual(solver, pose, x[i], y[i]);
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
    std::vector<Eigen::Vector2d> ori;
};

}
