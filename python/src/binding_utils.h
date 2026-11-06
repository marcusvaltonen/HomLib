#pragma once

#include <sstream>
#include <type_traits>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <Eigen/Dense>

#include <posedata.hpp>
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

    PySolver(HomLib::SolverId id, py::kwargs kwargs)
        : id(id), feature_type(HomLib::feature_type_for(id)), options(kwargs) {}
};

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

template <typename Solver>
Solver make_solver(const PySolver& py_solver) {
    Solver s{};
    HOMLIB_COPY_OPTION(s, py_solver.options, extra_check, bool)
    return s;
}
