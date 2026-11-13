#pragma once
#ifdef HOMLIB_HAS_MAGSAC

#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "magsac.h"
#include "model.h"
#include "model_score.h"
#include "samplers/uniform_sampler.h"
#include "samplers/prosac_sampler.h"
#include "samplers/importance_sampler.h"
#include "samplers/adaptive_reordering_sampler.h"

#include "posedata.hpp"
#include "gcransac_model_utils.h"
#include "gcransac_estimator.h"
#include "gcransac_solver.h"
#include "cv_data_utils.h"

namespace HomLib {

struct MAGSACStatistics {
    size_t num_iterations = 0;
    std::vector<size_t> inliers;
    double score = 0.0;
};

struct MAGSACSettings {
    double sigma_max = 10.0;
    double confidence = 0.99;
    size_t max_iterations = 5000;
    size_t min_iterations = 50;
    int partition_number = 5;
    bool use_magsac_plus_plus = true;
    GCSamplerType sampler = GCSamplerType::Uniform;
};

namespace detail {

template <class Estimator>
std::pair<PoseData, MAGSACStatistics> run_magsac_pipeline(
    Estimator& estimator,
    const cv::Mat& data,
    const std::vector<double>& probabilities,
    const MAGSACSettings& opts)
{
    using MAGSACType = MAGSAC<cv::Mat, Estimator>;
    using Version = typename MAGSACType::Version;
    using SamplerBase = gcransac::sampler::Sampler<cv::Mat, size_t>;

    MAGSACType magsac(opts.use_magsac_plus_plus
        ? Version::MAGSAC_PLUS_PLUS : Version::MAGSAC_ORIGINAL);
    magsac.setMaximumThreshold(opts.sigma_max);
    magsac.setIterationLimit(opts.max_iterations);
    magsac.setMinimumIterationNumber(opts.min_iterations);
    magsac.setPartitionNumber(opts.partition_number);
    magsac.setReferenceThreshold(opts.sigma_max / 2.0);

    const size_t sample_size = Estimator::sampleSize();
    if ((opts.sampler == GCSamplerType::ImportanceSampler ||
         opts.sampler == GCSamplerType::ARSampler) && probabilities.empty())
        throw std::runtime_error(
            "ImportanceSampler and ARSampler require probabilities.");
    std::unique_ptr<SamplerBase> sampler;
    switch (opts.sampler) {
        case GCSamplerType::PROSAC:
            sampler = std::make_unique<gcransac::sampler::ProsacSampler>(
                &data, sample_size);
            break;
        case GCSamplerType::ImportanceSampler:
            sampler = std::make_unique<gcransac::sampler::ImportanceSampler>(
                &data, probabilities, sample_size);
            break;
        case GCSamplerType::ARSampler:
            sampler = std::make_unique<gcransac::sampler::AdaptiveReorderingSampler>(
                &data, probabilities, sample_size);
            break;
        default:
            sampler = std::make_unique<gcransac::sampler::UniformSampler>(&data);
            break;
    }

    gcransac::Model model;
    int iterations = 0;
    ModelScore model_score;

    magsac.run(data, opts.confidence, estimator, *sampler,
               model, iterations, model_score);

    MAGSACStatistics stats;
    stats.num_iterations = static_cast<size_t>(iterations);
    stats.score = model_score.score;

    // Compute inliers from the best model
    PoseData pose;
    if (model.descriptor.size() > 0) {
        pose = gc_model_to_pose(model);
        const double sq_thresh = opts.sigma_max * opts.sigma_max;
        for (int i = 0; i < data.rows; ++i) {
            if (estimator.squaredResidual(data.row(i), model) < sq_thresh)
                stats.inliers.push_back(static_cast<size_t>(i));
        }
    }

    return {std::move(pose), std::move(stats)};
}

} // namespace detail

// ---------------------------------------------------------------------------
// Public API — run MAGSAC with HomLib solvers
// ---------------------------------------------------------------------------

template <class Solver>
std::pair<PoseData, MAGSACStatistics> run_magsac_point(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<double>& probabilities,
    const NonMinimalSettings& nm,
    const MAGSACSettings& settings,
    const NonMinimalEstimateFn& nm_solver_fn = {})
{
    constexpr size_t SS = Solver::minimal_sample_size_value;
    using MinSolver = GCPointMinimalSolver<Solver, SS>;
    using NonMinSolver = GCPointNonMinimalSolver<Solver, SS>;
    using Est = HomLibMAGSACEstimator<Solver, MinSolver, NonMinSolver>;

    auto data = detail::pack_point_cvmat(x, y);
    Est estimator(solver, nm, nm_solver_fn);
    auto result = detail::run_magsac_pipeline(estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

template <class Solver>
std::pair<PoseData, MAGSACStatistics> run_magsac_affine(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A,
    const std::vector<double>& probabilities,
    const NonMinimalSettings& nm,
    const MAGSACSettings& settings,
    const NonMinimalEstimateFn& nm_solver_fn = {})
{
    constexpr size_t SS = Solver::minimal_sample_size_value;
    using MinSolver = GCAffineMinimalSolver<Solver, SS>;
    using NonMinSolver = GCAffineNonMinimalSolver<Solver, SS>;
    using Est = HomLibMAGSACEstimator<Solver, MinSolver, NonMinSolver>;

    auto data = detail::pack_affine_cvmat(x, y, A);
    Est estimator(solver, nm, nm_solver_fn);
    auto result = detail::run_magsac_pipeline(estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

template <class Solver>
std::pair<PoseData, MAGSACStatistics> run_magsac_orientation(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori,
    const std::vector<double>& probabilities,
    const NonMinimalSettings& nm,
    const MAGSACSettings& settings,
    const NonMinimalEstimateFn& nm_solver_fn = {})
{
    constexpr size_t SS = Solver::minimal_sample_size_value;
    using MinSolver = GCOrientationMinimalSolver<Solver, SS>;
    using NonMinSolver = GCOrientationNonMinimalSolver<Solver, SS>;
    using Est = HomLibMAGSACEstimator<Solver, MinSolver, NonMinSolver>;

    auto data = detail::pack_orientation_cvmat(x, y, ori);
    Est estimator(solver, nm, nm_solver_fn);
    auto result = detail::run_magsac_pipeline(estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

} // namespace HomLib

#endif // HOMLIB_HAS_MAGSAC
