#include "binding_utils.h"
#include <magsac_runner.h>
#include <solver_dispatch.h>

using namespace pybind11::literals;

// --- Dispatch: MAGSAC ---

static std::pair<HomLib::PoseData, HomLib::MAGSACStatistics> dispatch_magsac_point(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const std::vector<double>& probabilities,
    HomLib::NonMinimalSettings nm,
    const py::object& non_minimal_solver_obj,
    const HomLib::MAGSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);

    if (!non_minimal_solver_obj.is_none())
        nm.solver_fn = resolve_non_minimal_fn(non_minimal_solver_obj.cast<PySolver>());

    return HomLib::dispatch_point_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_magsac_point(s, x, y, probabilities, nm, settings);
    });
}

static std::pair<HomLib::PoseData, HomLib::MAGSACStatistics> dispatch_magsac_affine(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic>& A_,
    const std::vector<double>& probabilities,
    const HomLib::MAGSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    std::vector<Eigen::Matrix2d> A;
    preprocess_affine(x_, y_, A_, &x, &y, &A);

    return HomLib::dispatch_affine_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_magsac_affine(s, x, y, A, probabilities, settings);
    });
}

static std::pair<HomLib::PoseData, HomLib::MAGSACStatistics> dispatch_magsac_orientation(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& ori_,
    const std::vector<double>& probabilities,
    const HomLib::MAGSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y, ori;
    preprocess_orientation(x_, y_, ori_, &x, &y, &ori);

    return HomLib::dispatch_orientation_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_magsac_orientation(s, x, y, ori, probabilities, settings);
    });
}

// --- pybind11 bindings ---

void bind_magsac(py::module_& m) {
    m.attr("_has_magsac") = true;

    py::class_<HomLib::MAGSACStatistics>(m, "MAGSACStatistics",
        "Statistics from a MAGSAC estimation run.")
        .def(py::init<>())
        .def_readwrite("num_iterations", &HomLib::MAGSACStatistics::num_iterations)
        .def_readwrite("inliers", &HomLib::MAGSACStatistics::inliers)
        .def_readwrite("score", &HomLib::MAGSACStatistics::score)
        .def("__repr__",
            [](const HomLib::MAGSACStatistics& s) {
                return "MAGSACStatistics("
                    "num_iterations=" + std::to_string(s.num_iterations) + ", "
                    "inliers=[list of length " + std::to_string(s.inliers.size()) + "], "
                    "score=" + std::to_string(s.score) + ")";
            });

    py::class_<HomLib::MAGSACSettings>(m, "MAGSACSettings",
        "Settings for MAGSAC estimation.")
        .def(py::init<>())
        .def_readwrite("sigma_max", &HomLib::MAGSACSettings::sigma_max)
        .def_readwrite("confidence", &HomLib::MAGSACSettings::confidence)
        .def_readwrite("max_iterations", &HomLib::MAGSACSettings::max_iterations)
        .def_readwrite("min_iterations", &HomLib::MAGSACSettings::min_iterations)
        .def_readwrite("partition_number", &HomLib::MAGSACSettings::partition_number)
        .def_readwrite("use_magsac_plus_plus", &HomLib::MAGSACSettings::use_magsac_plus_plus)
        .def_readwrite("sampler", &HomLib::MAGSACSettings::sampler)
        .def("__repr__",
            [](const HomLib::MAGSACSettings& s) {
                return "MAGSACSettings("
                    "sigma_max=" + std::to_string(s.sigma_max) + ", "
                    "confidence=" + std::to_string(s.confidence) + ", "
                    "max_iterations=" + std::to_string(s.max_iterations) + ")";
            });

    m.def("dispatch_magsac_point", &dispatch_magsac_point,
        "Run a point-based solver within the MAGSAC framework.",
        "solver"_a, "x"_a, "y"_a, "probabilities"_a = std::vector<double>{},
        "nm"_a = HomLib::NonMinimalSettings{},
        "non_minimal_solver"_a = py::none(), "settings"_a);
    m.def("dispatch_magsac_affine", &dispatch_magsac_affine,
        "Run an affine-covariant solver within the MAGSAC framework.",
        "solver"_a, "x"_a, "y"_a, "A"_a, "probabilities"_a = std::vector<double>{}, "settings"_a);
    m.def("dispatch_magsac_orientation", &dispatch_magsac_orientation,
        "Run an orientation-covariant solver within the MAGSAC framework.",
        "solver"_a, "x"_a, "y"_a, "ori"_a, "probabilities"_a = std::vector<double>{}, "settings"_a);
}
