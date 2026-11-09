#pragma once
#ifdef HOMLIB_HAS_GCRANSAC

#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "GCRANSAC.h"
#include "model.h"
#include "settings.h"
#include "statistics.h"
#include "scoring_function.h"
#include "backend_seed.h"
#include "samplers/uniform_sampler.h"
#include "neighborhood/grid_neighborhood_graph.h"
#include "preemption/preemption_empty.h"
#include "inlier_selectors/empty_inlier_selector.h"

#include "samplers/prosac_sampler.h"
#include "samplers/importance_sampler.h"
#include "samplers/adaptive_reordering_sampler.h"

#include "posedata.hpp"
#include "gcransac_model_utils.h"
#include "gcransac_estimator.h"
#include "gcransac_solver.h"
#include "cv_data_utils.h"

namespace HomLib {

struct GCRANSACStatistics {
    size_t num_iterations = 0;
    std::vector<size_t> inliers;
    double score = 0.0;
};

struct GCRANSACSettings {
    double threshold = 3.0;
    double confidence = 0.99;
    size_t max_iterations = 5000;
    size_t min_iterations = 50;
    double spatial_coherence_weight = 0.14;
    size_t cell_number = 4;
    bool do_local_optimization = true;
    bool do_graph_cut = true;
    bool do_final_least_squares = true;
    NonMinimalStrategy non_minimal_strategy = NonMinimalStrategy::SelectAndRefine;
    GCSamplerType sampler = GCSamplerType::Uniform;
    // The seed of the samplers (none: from std::random_device, a different run each time)
    std::optional<uint64_t> seed;
};

namespace detail {

template <size_t DimensionNumber, class Estimator>
std::pair<PoseData, GCRANSACStatistics> run_gcransac_pipeline(
    const Estimator& estimator,
    const cv::Mat& data,
    const std::vector<double>& probabilities,
    const GCRANSACSettings& opts)
{
    const backend_seed::Scope seed_scope(opts.seed);
    using NeighborhoodGraph = gcransac::neighborhood::GridNeighborhoodGraph<DimensionNumber>;
    using SamplerBase = gcransac::sampler::Sampler<cv::Mat, size_t>;

    // Build neighborhood graph.
    // GridNeighborhoodGraph computes cell indices as
    //   floor(coord / cell_size)  →  stored as size_t
    // so coordinates must be non-negative.  Shift each dimension
    // by its minimum so the grid sees values in [0, range].
    std::vector<double> cell_sizes(DimensionNumber);
    std::vector<double> lo_vals(DimensionNumber, 1e30);
    for (size_t i = 0; i < DimensionNumber; ++i) {
        double hi = -1e30;
        for (int r = 0; r < data.rows; ++r) {
            double v = data.at<double>(r, i);
            lo_vals[i] = std::min(lo_vals[i], v);
            hi = std::max(hi, v);
        }
        cell_sizes[i] = std::max(1.0, (hi - lo_vals[i]) / opts.cell_number);
    }

    cv::Mat shifted_data = data.clone();
    for (int r = 0; r < shifted_data.rows; ++r) {
        double* row = shifted_data.ptr<double>(r);
        for (size_t i = 0; i < DimensionNumber; ++i)
            row[i] -= lo_vals[i];
    }

    NeighborhoodGraph neighborhood(&shifted_data, cell_sizes, opts.cell_number);

    // Main sampler — type depends on settings + probabilities
    std::unique_ptr<SamplerBase> main_sampler;
    const size_t sample_size = Estimator::sampleSize();
    if ((opts.sampler == GCSamplerType::ImportanceSampler ||
         opts.sampler == GCSamplerType::ARSampler) && probabilities.empty())
        throw std::runtime_error(
            "ImportanceSampler and ARSampler require probabilities.");
    switch (opts.sampler) {
        case GCSamplerType::PROSAC:
            main_sampler = std::make_unique<gcransac::sampler::ProsacSampler>(
                &data, sample_size);
            break;
        case GCSamplerType::ImportanceSampler:
            main_sampler = std::make_unique<gcransac::sampler::ImportanceSampler>(
                &data, probabilities, sample_size);
            break;
        case GCSamplerType::ARSampler:
            main_sampler = std::make_unique<gcransac::sampler::AdaptiveReorderingSampler>(
                &data, probabilities, sample_size);
            break;
        default:
            main_sampler = std::make_unique<gcransac::sampler::UniformSampler>(&data);
            break;
    }
    gcransac::sampler::UniformSampler lo_sampler(&data);

    // Configure GCRANSAC
    using GC = gcransac::GCRANSAC<
        Estimator,
        NeighborhoodGraph,
        gcransac::MSACScoringFunction<Estimator>>;

    GC gcransac;
    gcransac.settings.threshold = opts.threshold;
    gcransac.settings.confidence = opts.confidence;
    gcransac.settings.max_iteration_number = opts.max_iterations;
    gcransac.settings.min_iteration_number = opts.min_iterations;
    gcransac.settings.spatial_coherence_weight = opts.spatial_coherence_weight;
    gcransac.settings.do_local_optimization = opts.do_local_optimization;
    gcransac.settings.do_graph_cut = opts.do_graph_cut;
    gcransac.settings.do_final_iterated_least_squares = opts.do_final_least_squares;

    gcransac::Model model;
    gcransac.run(data, estimator, main_sampler.get(), &lo_sampler,
                 &neighborhood, model);

    GCRANSACStatistics stats;
    const auto& gc_stats = gcransac.getRansacStatistics();
    stats.num_iterations = gc_stats.iteration_number;
    stats.inliers = gc_stats.inliers;
    stats.score = gc_stats.score;

    PoseData pose;
    if (!stats.inliers.empty())
        pose = gc_model_to_pose(model);
    return {std::move(pose), std::move(stats)};
}

} // namespace detail

// ---------------------------------------------------------------------------
// Public API — run GCRANSAC with HomLib solvers
// ---------------------------------------------------------------------------

template <class Solver>
std::pair<PoseData, GCRANSACStatistics> run_gcransac_point(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<double>& probabilities,
    const GCRANSACSettings& settings)
{
    constexpr size_t SS = Solver::minimal_sample_size_value;
    using MinSolver = GCPointMinimalSolver<Solver, SS>;
    using NonMinSolver = GCPointNonMinimalSolver<Solver, SS>;
    using Est = HomLibGCEstimator<Solver, MinSolver, NonMinSolver>;

    auto data = detail::pack_point_cvmat(x, y);
    Est estimator(solver, settings.non_minimal_strategy);
    auto result = detail::run_gcransac_pipeline<4>(estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

template <class Solver>
std::pair<PoseData, GCRANSACStatistics> run_gcransac_affine(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A,
    const std::vector<double>& probabilities,
    const GCRANSACSettings& settings)
{
    constexpr size_t SS = Solver::minimal_sample_size_value;
    using MinSolver = GCAffineMinimalSolver<Solver, SS>;
    using NonMinSolver = GCAffineNonMinimalSolver<Solver, SS>;
    using Est = HomLibGCEstimator<Solver, MinSolver, NonMinSolver>;

    auto data = detail::pack_affine_cvmat(x, y, A);
    Est estimator(solver, settings.non_minimal_strategy);
    auto result = detail::run_gcransac_pipeline<8>(estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

template <class Solver>
std::pair<PoseData, GCRANSACStatistics> run_gcransac_orientation(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori,
    const std::vector<double>& probabilities,
    const GCRANSACSettings& settings)
{
    constexpr size_t SS = Solver::minimal_sample_size_value;
    using MinSolver = GCOrientationMinimalSolver<Solver, SS>;
    using NonMinSolver = GCOrientationNonMinimalSolver<Solver, SS>;
    using Est = HomLibGCEstimator<Solver, MinSolver, NonMinSolver>;

    auto data = detail::pack_orientation_cvmat(x, y, ori);
    Est estimator(solver, settings.non_minimal_strategy);
    auto result = detail::run_gcransac_pipeline<6>(estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

} // namespace HomLib

#endif // HOMLIB_HAS_GCRANSAC
