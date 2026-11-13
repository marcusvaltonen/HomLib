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
#include <solver_utils.h>

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
        Eigen::Matrix2d tmp;   // the columns of A_ are row-major 2 x 2 matrices
        tmp << A_(0, i), A_(1, i),
               A_(2, i), A_(3, i);
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
    HomLib::DistortionCase distortion_case;
    bool supports_overdetermined;
    py::dict options;
    // The range on the distortion parameters set by Solver.set_k_range, in HomLib's normalized
    // coordinates (see DistortionRange); without it, the estimators' default, [-1, 0.1]
    bool has_k_range = false;
    // Whether the estimators normalize each sample (EstimatorBase::normalize_image_coord)
    bool normalize_image_coord = true;
    double k_min = -std::numeric_limits<double>::infinity();
    double k_max = std::numeric_limits<double>::infinity();

    PySolver(HomLib::SolverId id, py::kwargs kwargs)
        : id(id), feature_type(HomLib::feature_type_for(id)),
          distortion_case(HomLib::distortion_case_for(id)),
          supports_overdetermined(HomLib::supports_overdetermined_for(id)),
          options(kwargs) {}
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
HOMLIB_DEFINE_OPTION(use_qr)

// Detect flipped solver wrappers (they have a public inner_ member).
template <typename T, typename = void>
struct has_inner : std::false_type {};
template <typename T>
struct has_inner<T, std::void_t<decltype(std::declval<T&>().inner_)>> : std::true_type {};

template <typename Solver>
Solver make_solver(const PySolver& py_solver) {
    Solver s{};
    HOMLIB_COPY_OPTION(s, py_solver.options, extra_check, bool)
    HOMLIB_COPY_OPTION(s, py_solver.options, use_qr, bool)
    if constexpr (has_inner<Solver>::value) {
        HOMLIB_COPY_OPTION(s.inner_, py_solver.options, extra_check, bool)
        HOMLIB_COPY_OPTION(s.inner_, py_solver.options, use_qr, bool)
    }
    s.normalize_image_coord = py_solver.normalize_image_coord;
    if (py_solver.has_k_range) {
        s.k_range = user_range(py_solver);
    }
    return s;
}

inline HomLib::NonMinimalEstimateFn resolve_non_minimal_fn(const PySolver& py_solver) {
    switch (HomLib::feature_type_for(py_solver.id)) {
    case HomLib::FeatureType::Points:
        return HomLib::dispatch_point_solver(py_solver.id, [&py_solver](auto tag) -> HomLib::NonMinimalEstimateFn {
            using S = std::decay_t<decltype(tag)>;
            auto solver = make_solver<S>(py_solver);
            return [s = std::move(solver)](
                const std::vector<Eigen::Vector2d>& x,
                const std::vector<Eigen::Vector2d>& y,
                const std::vector<Eigen::Matrix2d>*,
                const std::vector<Eigen::Vector2d>*,
                std::vector<HomLib::PoseData>* poses) -> int {
                s.estimate(x, y, poses);
                return static_cast<int>(poses->size());
            };
        });
    case HomLib::FeatureType::Affine:
        return HomLib::dispatch_affine_solver(py_solver.id, [&py_solver](auto tag) -> HomLib::NonMinimalEstimateFn {
            using S = std::decay_t<decltype(tag)>;
            auto solver = make_solver<S>(py_solver);
            return [s = std::move(solver)](
                const std::vector<Eigen::Vector2d>& x,
                const std::vector<Eigen::Vector2d>& y,
                const std::vector<Eigen::Matrix2d>* A,
                const std::vector<Eigen::Vector2d>*,
                std::vector<HomLib::PoseData>* poses) -> int {
                if (!A) return 0;
                s.estimate(x, y, *A, poses);
                return static_cast<int>(poses->size());
            };
        });
    case HomLib::FeatureType::Orientation:
        return HomLib::dispatch_orientation_solver(py_solver.id, [&py_solver](auto tag) -> HomLib::NonMinimalEstimateFn {
            using S = std::decay_t<decltype(tag)>;
            auto solver = make_solver<S>(py_solver);
            return [s = std::move(solver)](
                const std::vector<Eigen::Vector2d>& x,
                const std::vector<Eigen::Vector2d>& y,
                const std::vector<Eigen::Matrix2d>*,
                const std::vector<Eigen::Vector2d>* ori,
                std::vector<HomLib::PoseData>* poses) -> int {
                if (!ori) return 0;
                s.estimate(x, y, *ori, poses);
                return static_cast<int>(poses->size());
            };
        });
    }
    return {};
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
