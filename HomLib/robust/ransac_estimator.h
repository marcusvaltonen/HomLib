#pragma once

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "solver_utils.h"

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
		if (nm_settings.sample_size > 0)
			return std::max(nm_settings.sample_size,
			                solver.minimal_sample_size() + 1);
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

	int NonMinimalSolver(const std::vector<int>& sample, HomLib::PoseData* pose) const {
		if (!use_non_minimal)
			return 0;

		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;
		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
		}

		auto est_fn = [](const Solver& s,
		                  const std::vector<Eigen::Vector2d>& xv,
		                  const std::vector<Eigen::Vector2d>& yv,
		                  const std::vector<size_t>& idx, size_t count,
		                  std::vector<HomLib::PoseData>* poses) {
			std::vector<Eigen::Vector2d> xs(count), ys(count);
			for (size_t i = 0; i < count; ++i) {
				xs[i] = xv[idx[i]];
				ys[i] = yv[idx[i]];
			}
			s.estimate(xs, ys, poses);
		};

		return detail::non_minimal_dispatch(
			solver, xx, yy, est_fn, *pose, nm_settings, detail::bind_nm(nm_solver_fn)) ? 1 : 0;
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
	NonMinimalSettings nm_settings;
	NonMinimalEstimateFn nm_solver_fn;
private:
	Solver solver;
	std::vector<Eigen::Vector2d> x;
	std::vector<Eigen::Vector2d> y;

};

}
