#pragma once

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "solver_utils.h"
namespace HomLib {

template<class Solver>
class AffineRansacEstimator {
public:
	AffineRansacEstimator(const std::vector<Eigen::Vector2d> &x_, const std::vector<Eigen::Vector2d> &y_, const std::vector<Eigen::Matrix2d> &A_, Solver est, const NonMinimalSettings& nm = {}) {
		x = x_;
		y = y_;
		A = A_;
		solver = est;
		nm_settings = nm;
	}

	inline int min_sample_size() const {
		return solver.minimal_sample_size();
	}
	inline int non_minimal_sample_size() const {
		if (nm_settings.sample_size > 0)
			return std::max(nm_settings.sample_size,
			                solver.minimal_sample_size() + 1);
		return solver.minimal_sample_size() * 4;
	}
	inline int num_data() const {
		return x.size();
	}

	int MinimalSolver(const std::vector<int>& sample,
		std::vector<HomLib::PoseData>* poses) const {
		
		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;
		std::vector<Eigen::Matrix2d> AA;

		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
			AA.push_back(A[sample[i]]);
		}
		solver.estimate(xx, yy, AA, poses);
		
		//std::cout << "H =\n" << (*poses)[0].homography << std::endl;
				
		return poses->size();
	}

	int NonMinimalSolver(const std::vector<int>& sample, HomLib::PoseData* pose) const {
		if (!use_non_minimal)
			return 0;

		std::vector<Eigen::Vector2d> xx;
		std::vector<Eigen::Vector2d> yy;
		std::vector<Eigen::Matrix2d> AA;
		for (size_t i = 0; i < sample.size(); i++) {
			xx.push_back(x[sample[i]]);
			yy.push_back(y[sample[i]]);
			AA.push_back(A[sample[i]]);
		}

		auto est_fn = [&AA](const Solver& s,
		                     const std::vector<Eigen::Vector2d>& xv,
		                     const std::vector<Eigen::Vector2d>& yv,
		                     const std::vector<size_t>& idx, size_t ss,
		                     std::vector<HomLib::PoseData>* poses) {
			std::vector<Eigen::Vector2d> xs(ss), ys(ss);
			std::vector<Eigen::Matrix2d> As(ss);
			for (size_t i = 0; i < ss; ++i) {
				xs[i] = xv[idx[i]]; ys[i] = yv[idx[i]]; As[i] = AA[idx[i]];
			}
			s.estimate(xs, ys, As, poses);
		};

		return detail::non_minimal_dispatch(
			solver, xx, yy, est_fn, *pose, nm_settings) ? 1 : 0;
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
private:
	Solver solver;
	std::vector<Eigen::Vector2d> x;
	std::vector<Eigen::Vector2d> y;
    std::vector<Eigen::Matrix2d> A;
};

}
