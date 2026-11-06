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
    py::dict options;
    // The range on the distortion parameters set by Solver.set_k_range, in HomLib's normalized
    // coordinates (see DistortionRange); without it, the estimators' default, [-1, 0.1]
    bool has_k_range = false;
    // Whether the estimators normalize each sample (EstimatorBase::normalize_image_coord)
    bool normalize_image_coord = true;
    double k_min = -std::numeric_limits<double>::infinity();
    double k_max = std::numeric_limits<double>::infinity();

    PySolver(HomLib::SolverId id, py::kwargs kwargs)
        : id(id), feature_type(HomLib::feature_type_for(id)), options(kwargs) {}
};

// The range set by Solver.set_k_range; set_k_range() without bounds means no check at all
inline HomLib::DistortionRange user_range(const PySolver& py_solver) {
    const double inf = std::numeric_limits<double>::infinity();
    if (py_solver.k_min == -inf && py_solver.k_max == inf) {
        return HomLib::DistortionRange();
    }
    return HomLib::DistortionRange(py_solver.k_min, py_solver.k_max);
}

// Stamp out a SFINAE member detector for a solver option.
#define HOMLIB_DEFINE_OPTION(name)                                                          \
    template <typename T, typename = void>                                                  \
    struct has_##name : std::false_type {};                                                  \
    template <typename T>                                                                   \
    struct has_##name<T, std::void_t<decltype(std::declval<T&>().name)>> : std::true_type {};

// Copy an option from PySolver.options into the C++ solver if the solver has that member.
#define HOMLIB_COPY_OPTION(solver, opts, name, type)                   \
    if constexpr (has_##name<std::decay_t<decltype(solver)>>::value) {  \
        if (opts.contains(#name))                                      \
            solver.name = opts[#name].cast<type>();                     \
    }

HOMLIB_DEFINE_OPTION(extra_check)

// Detect flipped solver wrappers (they have a public inner_ member).
template <typename T, typename = void>
struct has_inner : std::false_type {};
template <typename T>
struct has_inner<T, std::void_t<decltype(std::declval<T&>().inner_)>> : std::true_type {};

// The C++ solver of a PySolver: its options (extra_check), normalize_image_coord and, if set, its
// range on the distortion
template <typename Solver>
Solver make_solver(const PySolver& py_solver) {
    Solver s{};
    HOMLIB_COPY_OPTION(s, py_solver.options, extra_check, bool)
    if constexpr (has_inner<Solver>::value) {
        HOMLIB_COPY_OPTION(s.inner_, py_solver.options, extra_check, bool)
    }
    s.normalize_image_coord = py_solver.normalize_image_coord;
    if (py_solver.has_k_range) {
        s.k_range = user_range(py_solver);
    }
    return s;
}

// The range of estimate(), which calls the solver on the points as given (without the estimators'
// normalization, and without their default range): the range set by Solver.set_k_range, from the
// normalized coordinates to those of the points (k divided by f0^2, f0 = mean(|x|)/sqrt(2)), or
// as is without normalize_image_coord
inline HomLib::DistortionRange direct_range(const PySolver& py_solver, const std::vector<Eigen::Vector2d>& x) {
    const HomLib::DistortionRange range = py_solver.has_k_range ? user_range(py_solver) : HomLib::DistortionRange();
    if (!range.enabled || !py_solver.normalize_image_coord || x.empty()) {
        return range;
    }
    double f0 = 0.0;
    for (const Eigen::Vector2d& p : x) {
        f0 += p.norm();
    }
    f0 /= x.size() * std::sqrt(2.0);
    return HomLib::DistortionRange(py_solver.k_min / (f0 * f0), py_solver.k_max / (f0 * f0));
}
