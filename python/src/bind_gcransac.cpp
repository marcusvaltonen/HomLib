#include "binding_utils.h"
#include <gcransac_runner.h>
#include <solver_dispatch.h>

using namespace pybind11::literals;

// --- Dispatch: GCRANSAC ---

static std::pair<HomLib::PoseData, HomLib::GCRANSACStatistics> dispatch_gcransac_point(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const std::vector<double>& probabilities,
    HomLib::NonMinimalSettings nm,
    const py::object& non_minimal_solver_obj,
    const HomLib::GCRANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);

    HomLib::NonMinimalEstimateFn nm_solver_fn;
    if (!non_minimal_solver_obj.is_none())
        nm_solver_fn = resolve_non_minimal_fn(non_minimal_solver_obj.cast<PySolver>());

    return HomLib::dispatch_point_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_gcransac_point(s, x, y, probabilities, nm, settings, nm_solver_fn);
    });
}

static std::pair<HomLib::PoseData, HomLib::GCRANSACStatistics> dispatch_gcransac_affine(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic>& A_,
    const std::vector<double>& probabilities,
    HomLib::NonMinimalSettings nm,
    const py::object& non_minimal_solver_obj,
    const HomLib::GCRANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    std::vector<Eigen::Matrix2d> A;
    preprocess_affine(x_, y_, A_, &x, &y, &A);

    HomLib::NonMinimalEstimateFn nm_solver_fn;
    if (!non_minimal_solver_obj.is_none())
        nm_solver_fn = resolve_non_minimal_fn(non_minimal_solver_obj.cast<PySolver>());

    return HomLib::dispatch_affine_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_gcransac_affine(s, x, y, A, probabilities, nm, settings, nm_solver_fn);
    });
}

static std::pair<HomLib::PoseData, HomLib::GCRANSACStatistics> dispatch_gcransac_orientation(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& ori_,
    const std::vector<double>& probabilities,
    HomLib::NonMinimalSettings nm,
    const py::object& non_minimal_solver_obj,
    const HomLib::GCRANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y, ori;
    preprocess_orientation(x_, y_, ori_, &x, &y, &ori);

    HomLib::NonMinimalEstimateFn nm_solver_fn;
    if (!non_minimal_solver_obj.is_none())
        nm_solver_fn = resolve_non_minimal_fn(non_minimal_solver_obj.cast<PySolver>());

    return HomLib::dispatch_orientation_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_gcransac_orientation(s, x, y, ori, probabilities, nm, settings, nm_solver_fn);
    });
}

// --- pybind11 bindings ---

void bind_gcransac(py::module_& m) {
    m.attr("_has_gcransac") = true;

    py::class_<HomLib::GCRANSACStatistics>(m, "GCRANSACStatistics",
        "Statistics from a GC-RANSAC estimation run.")
        .def(py::init<>())
        .def_readwrite("num_iterations", &HomLib::GCRANSACStatistics::num_iterations)
        .def_readwrite("inliers", &HomLib::GCRANSACStatistics::inliers)
        .def_readwrite("score", &HomLib::GCRANSACStatistics::score)
        .def("__repr__",
            [](const HomLib::GCRANSACStatistics& s) {
                return "GCRANSACStatistics("
                    "num_iterations=" + std::to_string(s.num_iterations) + ", "
                    "inliers=[list of length " + std::to_string(s.inliers.size()) + "], "
                    "score=" + std::to_string(s.score) + ")";
            });

    py::class_<HomLib::GCRANSACSettings>(m, "GCRANSACSettings",
        "Settings for GC-RANSAC estimation.")
        .def(py::init<>())
        .def_readwrite("threshold", &HomLib::GCRANSACSettings::threshold)
        .def_readwrite("confidence", &HomLib::GCRANSACSettings::confidence)
        .def_readwrite("max_iterations", &HomLib::GCRANSACSettings::max_iterations)
        .def_readwrite("min_iterations", &HomLib::GCRANSACSettings::min_iterations)
        .def_readwrite("spatial_coherence_weight", &HomLib::GCRANSACSettings::spatial_coherence_weight)
        .def_readwrite("cell_number", &HomLib::GCRANSACSettings::cell_number)
        .def_readwrite("do_local_optimization", &HomLib::GCRANSACSettings::do_local_optimization)
        .def_readwrite("do_graph_cut", &HomLib::GCRANSACSettings::do_graph_cut)
        .def_readwrite("do_final_least_squares", &HomLib::GCRANSACSettings::do_final_least_squares)
        .def_readwrite("sampler", &HomLib::GCRANSACSettings::sampler)
        .def("__repr__",
            [](const HomLib::GCRANSACSettings& s) {
                return "GCRANSACSettings("
                    "threshold=" + std::to_string(s.threshold) + ", "
                    "confidence=" + std::to_string(s.confidence) + ", "
                    "max_iterations=" + std::to_string(s.max_iterations) + ")";
            });

    m.def("dispatch_gcransac_point", &dispatch_gcransac_point,
        "Run a point-based solver within the GC-RANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "probabilities"_a = std::vector<double>{},
        "nm"_a = HomLib::NonMinimalSettings{},
        "non_minimal_solver"_a = py::none(), "settings"_a);
    m.def("dispatch_gcransac_affine", &dispatch_gcransac_affine,
        "Run an affine-covariant solver within the GC-RANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "A"_a, "probabilities"_a = std::vector<double>{},
        "nm"_a = HomLib::NonMinimalSettings{},
        "non_minimal_solver"_a = py::none(), "settings"_a);
    m.def("dispatch_gcransac_orientation", &dispatch_gcransac_orientation,
        "Run an orientation-covariant solver within the GC-RANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "ori"_a, "probabilities"_a = std::vector<double>{},
        "nm"_a = HomLib::NonMinimalSettings{},
        "non_minimal_solver"_a = py::none(), "settings"_a);
}
