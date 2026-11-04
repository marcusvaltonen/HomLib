#pragma once

#include <cmath>
#include <limits>
#include <optional>

#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <Eigen/Dense>

#include <posedata.hpp>
#include <distortion_range.h>
#include <solver_id.h>
#include <solver_dispatch.h>

namespace py = pybind11;

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

inline void preprocess(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    std::vector<Eigen::Vector2d> *x,
    std::vector<Eigen::Vector2d> *y
) {
    if (x_.cols() != y_.cols()) {
        throw std::invalid_argument("x and y should be of equal size.");
    }
    for (Eigen::Index i=0; i < x_.cols(); i++) {
        x->push_back(x_.col(i));
        y->push_back(y_.col(i));
    }
}

inline void preprocess_affine(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic> &A_,
    std::vector<Eigen::Vector2d> *x,
    std::vector<Eigen::Vector2d> *y,
    std::vector<Eigen::Matrix2d> *A
) {
    if (x_.cols() != y_.cols()) {
        throw std::invalid_argument("x and y should be of equal size.");
    }
    if (x_.cols() != A_.cols()) {
        throw std::invalid_argument("x and A should be of equal size.");
    }
    for (size_t i=0; i < x_.cols(); i++) {
        x->push_back(x_.col(i));
        y->push_back(y_.col(i));
        Eigen::Matrix2d tmp;
        tmp << A_.col(i);
        tmp.transposeInPlace();
        A->push_back(tmp);
    }
}

inline void preprocess_orientation(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &ori_,
    std::vector<Eigen::Vector2d> *x,
    std::vector<Eigen::Vector2d> *y,
    std::vector<Eigen::Vector2d> *ori
) {
    if (x_.cols() != y_.cols()) {
        throw std::invalid_argument("x and y should be of equal size.");
    }
    if (x_.cols() != ori_.cols()) {
        throw std::invalid_argument("x and ori should be of equal size.");
    }
    for (size_t i=0; i < x_.cols(); i++) {
        x->push_back(x_.col(i));
        y->push_back(y_.col(i));
        ori->push_back(ori_.col(i));
    }
}

struct PySolver {
    HomLib::SolverId id;
    HomLib::FeatureType feature_type;
    bool extra_check;
    // Optional range on the distortion parameters of the solutions (Solver.set_k_range), in the
    // units of the points given to estimate(); without it, the solver's default range applies
    // (PoseEstimator::use_k_range: on for the two-sided Wadenbäck solver, off otherwise)
    bool has_k_range = false;
    double k_min = -std::numeric_limits<double>::infinity();
    double k_max = std::numeric_limits<double>::infinity();

    PySolver(HomLib::SolverId id, bool extra_check = false)
        : id(id), feature_type(HomLib::feature_type_for(id)), extra_check(extra_check) {}
};

// The C++ solver of a PySolver: its extra_check and, if set, its range on the distortion
template <typename Solver>
Solver make_solver(const PySolver& py_solver) {
    Solver s = [&]() {
        if constexpr (std::is_constructible_v<Solver, bool>) {
            return Solver(py_solver.extra_check);
        } else {
            return Solver();
        }
    }();
    if (py_solver.has_k_range) {
        // without bounds (set_k_range()): no check at all, not even the solver's default
        s.use_k_range = std::isfinite(py_solver.k_min) || std::isfinite(py_solver.k_max);
        s.k_range_in_input_units = true;
        s.k_min = py_solver.k_min;
        s.k_max = py_solver.k_max;
    }
    return s;
}

// The range that solve() applies when estimate() calls the solver directly (without the
// estimator's normalization): the solver's range, applied to the points as given. For the
// default range of the two-sided Wadenbäck solver, [-1, 0.1], this is the check its solver
// function made before the range moved to the estimators.
template <typename Solver>
HomLib::DistortionRange direct_range(const Solver& s) {
    return s.use_k_range ? HomLib::DistortionRange(s.k_min, s.k_max) : HomLib::DistortionRange();
}
