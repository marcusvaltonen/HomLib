#include "binding_utils.h"
#include <superansac_runner.h>

using namespace pybind11::literals;

// --- Dispatch: SupeRANSAC ---

static std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_point(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);

    return HomLib::dispatch_point_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_superansac_point(s, x, y, settings);
    });
}

static std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_affine(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic>& A_,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    std::vector<Eigen::Matrix2d> A;
    preprocess_affine(x_, y_, A_, &x, &y, &A);

    return HomLib::dispatch_affine_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_superansac_affine(s, x, y, A, settings);
    });
}

static std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_orientation(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& ori_,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y, ori;
    preprocess_orientation(x_, y_, ori_, &x, &y, &ori);

    return HomLib::dispatch_orientation_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_superansac_orientation(s, x, y, ori, settings);
    });
}

// --- pybind11 bindings ---

void bind_superansac(py::module_& m) {
    m.attr("_has_superansac") = true;

    py::class_<HomLib::SuperansacStatistics>(m, "SuperansacStatistics",
        "Statistics from a SupeRANSAC estimation run.")
        .def(py::init<>())
        .def_readwrite("num_iterations", &HomLib::SuperansacStatistics::num_iterations,
            "Total number of RANSAC iterations.")
        .def_readwrite("inliers", &HomLib::SuperansacStatistics::inliers,
            "Indices of inlier correspondences.")
        .def_readwrite("best_score", &HomLib::SuperansacStatistics::best_score,
            "Score of the best model.")
        .def("__repr__",
            [](const HomLib::SuperansacStatistics& s) {
                return "SuperansacStatistics("
                    "num_iterations=" + std::to_string(s.num_iterations) + ", "
                    "inliers=[list of length " + std::to_string(s.inliers.size()) + "], "
                    "best_score=" + std::to_string(s.best_score) +
                    ")";
            });

    py::class_<superansac::RANSACSettings>(m, "SuperansacSettings",
        "Settings for SupeRANSAC estimation.")
        .def(py::init<>())
        .def_readwrite("inlier_threshold", &superansac::RANSACSettings::inlierThreshold,
            "Inlier threshold.")
        .def_readwrite("confidence", &superansac::RANSACSettings::confidence,
            "Desired confidence level.")
        .def_readwrite("min_iterations", &superansac::RANSACSettings::minIterations,
            "Minimum number of iterations.")
        .def_readwrite("max_iterations", &superansac::RANSACSettings::maxIterations,
            "Maximum number of iterations.")
        .def_readwrite("top_k", &superansac::RANSACSettings::topKForLocalOptimization,
            "Number of top models for local optimization.")
        .def_readwrite("use_sprt", &superansac::RANSACSettings::useSprt,
            "Whether to use SPRT test.")
        .def_readwrite("sampler", &superansac::RANSACSettings::sampler,
            "Sampler type (default: PROSAC). Use SuperansacSamplerType enum.")
        .def_readwrite("scoring", &superansac::RANSACSettings::scoring,
            "Scoring type (default: MAGSAC). Use SuperansacScoringType enum.")
        .def_readwrite("local_optimization", &superansac::RANSACSettings::localOptimization,
            "Local optimization type. Use SuperansacLocalOptType enum.")
        .def_readwrite("final_optimization", &superansac::RANSACSettings::finalOptimization,
            "Final optimization type. Use SuperansacLocalOptType enum.")
        .def_readwrite("neighborhood", &superansac::RANSACSettings::neighborhood,
            "Neighborhood type (Grid, FLANN_KNN, FLANN_Radius). Use SuperansacNeighborhoodType enum.")
        .def_readwrite("neighborhood_settings", &superansac::RANSACSettings::neighborhoodSettings,
            "Neighborhood settings (grid density, size, etc.).")
        .def("__repr__",
            [](const superansac::RANSACSettings& s) {
                return "SuperansacSettings("
                    "inlier_threshold=" + std::to_string(s.inlierThreshold) + ", "
                    "confidence=" + std::to_string(s.confidence) + ", "
                    "min_iterations=" + std::to_string(s.minIterations) + ", "
                    "max_iterations=" + std::to_string(s.maxIterations) +
                    ")";
            });

    py::enum_<superansac::scoring::ScoringType>(m, "SuperansacScoringType")
        .value("RANSAC", superansac::scoring::ScoringType::RANSAC)
        .value("MSAC", superansac::scoring::ScoringType::MSAC)
        .value("MAGSAC", superansac::scoring::ScoringType::MAGSAC)
        .value("MINPRAN", superansac::scoring::ScoringType::MINPRAN)
        .value("ACRANSAC", superansac::scoring::ScoringType::ACRANSAC)
        .value("GAU", superansac::scoring::ScoringType::GAU)
        .value("ML", superansac::scoring::ScoringType::ML)
        .value("GRID", superansac::scoring::ScoringType::GRID);

    py::enum_<superansac::samplers::SamplerType>(m, "SuperansacSamplerType")
        .value("Uniform", superansac::samplers::SamplerType::Uniform)
        .value("PROSAC", superansac::samplers::SamplerType::PROSAC)
        .value("NAPSAC", superansac::samplers::SamplerType::NAPSAC)
        .value("ProgressiveNAPSAC", superansac::samplers::SamplerType::ProgressiveNAPSAC)
        .value("ImportanceSampler", superansac::samplers::SamplerType::ImportanceSampler)
        .value("ARSampler", superansac::samplers::SamplerType::ARSampler);

    py::class_<superansac::NeighborhoodSettings>(m, "SuperansacNeighborhoodSettings",
        "Settings for the spatial neighborhood graph.")
        .def(py::init<>())
        .def_readwrite("neighborhood_size", &superansac::NeighborhoodSettings::neighborhoodSize,
            "Neighborhood size (radius for FLANN, ignored for Grid).")
        .def_readwrite("neighborhood_grid_density", &superansac::NeighborhoodSettings::neighborhoodGridDensity,
            "Number of cells along each axis for Grid neighborhood.")
        .def_readwrite("nearest_neighbor_number", &superansac::NeighborhoodSettings::nearestNeighborNumber,
            "Number of nearest neighbors for FLANN_KNN.");

    py::enum_<superansac::neighborhood::NeighborhoodType>(m, "SuperansacNeighborhoodType")
        .value("Grid", superansac::neighborhood::NeighborhoodType::Grid)
        .value("FLANN_KNN", superansac::neighborhood::NeighborhoodType::FLANN_KNN)
        .value("FLANN_Radius", superansac::neighborhood::NeighborhoodType::FLANN_Radius);

    py::enum_<superansac::local_optimization::LocalOptimizationType>(m, "SuperansacLocalOptType")
        .value("None_", superansac::local_optimization::LocalOptimizationType::None)
        .value("LSQ", superansac::local_optimization::LocalOptimizationType::LSQ)
        .value("IRLS", superansac::local_optimization::LocalOptimizationType::IRLS)
        .value("NestedRANSAC", superansac::local_optimization::LocalOptimizationType::NestedRANSAC)
        .value("GCRANSAC", superansac::local_optimization::LocalOptimizationType::GCRANSAC)
        .value("IteratedLMEDS", superansac::local_optimization::LocalOptimizationType::IteratedLMEDS)
        .value("CrossValidation", superansac::local_optimization::LocalOptimizationType::CrossValidation);

    m.def("dispatch_superansac_point", &dispatch_superansac_point,
        "Run a point-based solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "settings"_a);

    m.def("dispatch_superansac_affine", &dispatch_superansac_affine,
        "Run an affine-covariant solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "A"_a, "settings"_a);

    m.def("dispatch_superansac_orientation", &dispatch_superansac_orientation,
        "Run an orientation-covariant solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "ori"_a, "settings"_a);
}
