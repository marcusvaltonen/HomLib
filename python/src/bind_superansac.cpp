#include "binding_utils.h"
#include <superansac_runner.h>
#include <solver_dispatch.h>

using namespace pybind11::literals;

// --- Dispatch: SupeRANSAC ---

static std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_point(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const std::vector<double>& probabilities,
    HomLib::NonMinimalSettings nm,
    const py::object& non_minimal_solver_obj,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);

    HomLib::NonMinimalEstimateFn nm_solver_fn;
    if (!non_minimal_solver_obj.is_none())
        nm_solver_fn = resolve_non_minimal_fn(non_minimal_solver_obj.cast<PySolver>());

    return HomLib::dispatch_point_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_superansac_point(s, x, y, probabilities, nm, settings, nm_solver_fn);
    });
}

static std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_affine(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic>& A_,
    const std::vector<double>& probabilities,
    HomLib::NonMinimalSettings nm,
    const py::object& non_minimal_solver_obj,
    const superansac::RANSACSettings& settings)
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
        return HomLib::run_superansac_affine(s, x, y, A, probabilities, nm, settings, nm_solver_fn);
    });
}

static std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_orientation(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& ori_,
    const std::vector<double>& probabilities,
    HomLib::NonMinimalSettings nm,
    const py::object& non_minimal_solver_obj,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y, ori;
    preprocess_orientation(x_, y_, ori_, &x, &y, &ori);

    HomLib::NonMinimalEstimateFn nm_solver_fn;
    if (!non_minimal_solver_obj.is_none())
        nm_solver_fn = resolve_non_minimal_fn(non_minimal_solver_obj.cast<PySolver>());

    return HomLib::dispatch_orientation_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver);
        return HomLib::run_superansac_orientation(s, x, y, ori, probabilities, nm, settings, nm_solver_fn);
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

    // --- SupeRANSAC sub-settings structs ---
    py::class_<superansac::ARSamplerSettings>(m, "SuperansacARSamplerSettings")
        .def(py::init<>())
        .def_readwrite("estimator_variance", &superansac::ARSamplerSettings::estimatorVariance)
        .def_readwrite("randomness", &superansac::ARSamplerSettings::randomness);

    py::class_<superansac::LocalOptimizationSettings>(m, "SuperansacLOSettings")
        .def(py::init<>())
        .def_readwrite("max_iterations", &superansac::LocalOptimizationSettings::maxIterations)
        .def_readwrite("graph_cut_number", &superansac::LocalOptimizationSettings::graphCutNumber)
        .def_readwrite("sample_size_multiplier", &superansac::LocalOptimizationSettings::sampleSizeMultiplier)
        .def_readwrite("spatial_coherence_weight", &superansac::LocalOptimizationSettings::spatialCoherenceWeight);

    py::class_<superansac::DegensacSettings>(m, "SuperansacDegensacSettings")
        .def(py::init<>())
        .def_readwrite("enabled", &superansac::DegensacSettings::enabled)
        .def_readwrite("homography_threshold", &superansac::DegensacSettings::homographyThreshold)
        .def_readwrite("inner_ransac_iterations", &superansac::DegensacSettings::innerRansacIterations);

    py::class_<superansac::NeighborhoodSettings>(m, "SuperansacNeighborhoodSettings")
        .def(py::init<>())
        .def_readwrite("neighborhood_size", &superansac::NeighborhoodSettings::neighborhoodSize)
        .def_readwrite("neighborhood_grid_density", &superansac::NeighborhoodSettings::neighborhoodGridDensity)
        .def_readwrite("nearest_neighbor_number", &superansac::NeighborhoodSettings::nearestNeighborNumber);

    // --- SupeRANSAC main settings ---
    py::class_<superansac::RANSACSettings>(m, "SuperansacSettings",
        "Settings for SupeRANSAC estimation.")
        .def(py::init<>())
        .def_readwrite("inlier_threshold", &superansac::RANSACSettings::inlierThreshold)
        .def_readwrite("confidence", &superansac::RANSACSettings::confidence)
        .def_readwrite("min_iterations", &superansac::RANSACSettings::minIterations)
        .def_readwrite("max_iterations", &superansac::RANSACSettings::maxIterations)
        .def_readwrite("top_k", &superansac::RANSACSettings::topKForLocalOptimization)
        .def_readwrite("use_sprt", &superansac::RANSACSettings::useSprt)
        .def_readwrite("sampler", &superansac::RANSACSettings::sampler)
        .def_readwrite("scoring", &superansac::RANSACSettings::scoring)
        .def_readwrite("local_optimization", &superansac::RANSACSettings::localOptimization)
        .def_readwrite("final_optimization", &superansac::RANSACSettings::finalOptimization)
        .def_readwrite("neighborhood", &superansac::RANSACSettings::neighborhood)
        .def_readwrite("termination_criterion", &superansac::RANSACSettings::terminationCriterion)
        .def_readwrite("inlier_selector", &superansac::RANSACSettings::inlierSelector)
        // Bool flags
        .def_readwrite("local_optimization_inside_loop", &superansac::RANSACSettings::localOptimizationInsideTheLoop)
        .def_readwrite("fundamental_orientation_check", &superansac::RANSACSettings::fundamentalOrientationCheck)
        .def_readwrite("fundamental_symmetric_validation", &superansac::RANSACSettings::fundamentalSymmetricValidation)
        .def_readwrite("final_polish_top_k", &superansac::RANSACSettings::finalPolishTopK)
        .def_readwrite("final_refinement_weighted", &superansac::RANSACSettings::finalRefinementWeighted)
        .def_readwrite("homography_bundle_refinement", &superansac::RANSACSettings::homographyBundleRefinement)
        // Scalar settings
        .def_readwrite("final_refinement_rounds", &superansac::RANSACSettings::finalRefinementRounds)
        .def_readwrite("final_refinement_schedule_start", &superansac::RANSACSettings::finalRefinementScheduleStart)
        .def_readwrite("final_refinement_schedule_end", &superansac::RANSACSettings::finalRefinementScheduleEnd)
        .def_readwrite("final_refinement_schedule", &superansac::RANSACSettings::finalRefinementSchedule)
        .def_readwrite("essential_cheirality_point_number", &superansac::RANSACSettings::essentialCheiralityPointNumber)
        .def_readwrite("homography_bundle_max_iterations", &superansac::RANSACSettings::homographyBundleMaxIterations)
        // Sub-settings structs
        .def_readwrite("ar_sampler_settings", &superansac::RANSACSettings::arSamplerSettings)
        .def_readwrite("lo_settings", &superansac::RANSACSettings::localOptimizationSettings)
        .def_readwrite("final_opt_settings", &superansac::RANSACSettings::finalOptimizationSettings)
        .def_readwrite("neighborhood_settings", &superansac::RANSACSettings::neighborhoodSettings)
        .def_readwrite("degensac_settings", &superansac::RANSACSettings::degensacSettings)
        .def("__repr__",
            [](const superansac::RANSACSettings& s) {
                return "SuperansacSettings("
                    "inlier_threshold=" + std::to_string(s.inlierThreshold) + ", "
                    "confidence=" + std::to_string(s.confidence) + ", "
                    "min_iterations=" + std::to_string(s.minIterations) + ", "
                    "max_iterations=" + std::to_string(s.maxIterations) +
                    ")";
            });

    // --- SupeRANSAC enums ---
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

    py::enum_<superansac::termination::TerminationType>(m, "SuperansacTerminationType")
        .value("RANSAC", superansac::termination::TerminationType::RANSAC);

    py::enum_<superansac::inlier_selector::InlierSelectorType>(m, "SuperansacInlierSelectorType")
        .value("None_", superansac::inlier_selector::InlierSelectorType::None)
        .value("SpacePartitioningRANSAC", superansac::inlier_selector::InlierSelectorType::SpacePartitioningRANSAC);

    m.def("dispatch_superansac_point", &dispatch_superansac_point,
        "Run a point-based solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "probabilities"_a = std::vector<double>{},
        "nm"_a = HomLib::NonMinimalSettings{},
        "non_minimal_solver"_a = py::none(), "settings"_a);

    m.def("dispatch_superansac_affine", &dispatch_superansac_affine,
        "Run an affine-covariant solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "A"_a, "probabilities"_a = std::vector<double>{},
        "nm"_a = HomLib::NonMinimalSettings{},
        "non_minimal_solver"_a = py::none(), "settings"_a);

    m.def("dispatch_superansac_orientation", &dispatch_superansac_orientation,
        "Run an orientation-covariant solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "ori"_a, "probabilities"_a = std::vector<double>{},
        "nm"_a = HomLib::NonMinimalSettings{},
        "non_minimal_solver"_a = py::none(), "settings"_a);
}
