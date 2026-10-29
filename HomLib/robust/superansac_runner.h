#pragma once
#ifdef HOMLIB_HAS_SUPERANSAC

#include <utility>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include "superansac.h"
#include "settings.h"
#include "samplers/types.h"
#include "samplers/prosac_sampler.h"
#include "scoring/types.h"
#include "local_optimization/types.h"
#include "local_optimization/irls_optimizer.h"
#include "local_optimization/least_squares_optimizer.h"
#include "local_optimization/nested_ransac_optimizer.h"
#include "termination/types.h"
#include "termination/ransac_criterion.h"
#include "scoring/magsac_scoring.h"
#include "scoring/magsac_sprt_scoring.h"
#include "neighborhood/types.h"
#include "inlier_selectors/types.h"

#include "posedata.hpp"
#include "superansac_model_utils.h"
#include "superansac_estimator.h"

namespace HomLib {

struct SuperansacStatistics {
    size_t num_iterations = 0;
    std::vector<size_t> inliers;
    double best_score = 0.0;
};

namespace detail {

// Pack point correspondences into DataMatrix.
inline DataMatrix pack_point_data(
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y)
{
    const size_t n = x.size();
    DataMatrix data(n, 4);
    for (size_t i = 0; i < n; ++i) {
        data(i, 0) = x[i](0);
        data(i, 1) = x[i](1);
        data(i, 2) = y[i](0);
        data(i, 3) = y[i](1);
    }
    return data;
}

inline DataMatrix pack_affine_data(
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A)
{
    const size_t n = x.size();
    DataMatrix data(n, 8);
    for (size_t i = 0; i < n; ++i) {
        data(i, 0) = x[i](0);
        data(i, 1) = x[i](1);
        data(i, 2) = y[i](0);
        data(i, 3) = y[i](1);
        data(i, 4) = A[i](0, 0);
        data(i, 5) = A[i](0, 1);
        data(i, 6) = A[i](1, 0);
        data(i, 7) = A[i](1, 1);
    }
    return data;
}

inline DataMatrix pack_orientation_data(
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori)
{
    const size_t n = x.size();
    DataMatrix data(n, 6);
    for (size_t i = 0; i < n; ++i) {
        data(i, 0) = x[i](0);
        data(i, 1) = x[i](1);
        data(i, 2) = y[i](0);
        data(i, 3) = y[i](1);
        data(i, 4) = ori[i](0);
        data(i, 5) = ori[i](1);
    }
    return data;
}

inline void initialize_local_optimizer(
    std::unique_ptr<superansac::local_optimization::LocalOptimizer>& opt,
    superansac::local_optimization::LocalOptimizationType type,
    const superansac::LocalOptimizationSettings& loSettings,
    bool isFinal)
{
    using namespace superansac::local_optimization;
    switch (type) {
        case LocalOptimizationType::IRLS: {
            auto* irls = dynamic_cast<IRLSOptimizer*>(opt.get());
            irls->setMaxIterations(loSettings.maxIterations);
            irls->setUseInliers(isFinal || true);
            break;
        }
        case LocalOptimizationType::LSQ: {
            auto* lsq = dynamic_cast<LeastSquaresOptimizer*>(opt.get());
            lsq->setUseInliers(isFinal || true);
            break;
        }
        case LocalOptimizationType::NestedRANSAC: {
            auto* nested = dynamic_cast<NestedRANSACOptimizer*>(opt.get());
            nested->setMaxIterations(loSettings.maxIterations);
            nested->setSampleSizeMultiplier(loSettings.sampleSizeMultiplier);
            break;
        }
        default:
            break;
    }
}

// Assemble and run the SupeRANSAC pipeline with a given estimator and data.
// DimensionNumber is the number of columns in the DataMatrix (4, 6, or 8).
template <size_t DimensionNumber>
std::pair<PoseData, SuperansacStatistics> run_pipeline(
    superansac::estimator::Estimator* estimator,
    const DataMatrix& data,
    const superansac::RANSACSettings& settings)
{
    auto sampler = superansac::samplers::createSampler<DimensionNumber>(settings.sampler);
    auto scoring = superansac::scoring::createScoring<DimensionNumber>(settings.scoring, settings.useSprt);
    auto termination = superansac::termination::createTerminationCriterion(settings.terminationCriterion);

    // Sampler initialization (required before run)
    if (settings.sampler == superansac::samplers::SamplerType::PROSAC)
        dynamic_cast<superansac::samplers::PROSACSampler*>(sampler.get())
            ->setSampleSize(estimator->sampleSize());
    sampler->initialize(data);

    // Scoring initialization
    scoring->setThreshold(settings.inlierThreshold);
    if (settings.scoring == superansac::scoring::ScoringType::MAGSAC) {
        if (settings.useSprt)
            dynamic_cast<superansac::scoring::MAGSACSPRTScoring*>(scoring.get())->initialize(estimator);
        else
            dynamic_cast<superansac::scoring::MAGSACScoring*>(scoring.get())->initialize(estimator);
    }

    // Termination criterion initialization
    if (settings.terminationCriterion == superansac::termination::TerminationType::RANSAC)
        dynamic_cast<superansac::termination::RANSACCriterion*>(termination.get())
            ->setConfidence(settings.confidence);

    // Local and final optimizer initialization
    std::unique_ptr<superansac::local_optimization::LocalOptimizer> localOpt, finalOpt;
    if (settings.localOptimization != superansac::local_optimization::LocalOptimizationType::None) {
        localOpt = superansac::local_optimization::createLocalOptimizer(settings.localOptimization);
        initialize_local_optimizer(localOpt, settings.localOptimization, settings.localOptimizationSettings, false);
    }
    if (settings.finalOptimization != superansac::local_optimization::LocalOptimizationType::None) {
        finalOpt = superansac::local_optimization::createLocalOptimizer(settings.finalOptimization);
        initialize_local_optimizer(finalOpt, settings.finalOptimization, settings.finalOptimizationSettings, true);
    }

    superansac::SupeRansac ransac;
    ransac.setEstimator(estimator);
    ransac.setSampler(sampler.get());
    ransac.setScoring(scoring.get());
    ransac.setTerminationCriterion(termination.get());
    if (localOpt) ransac.setLocalOptimizer(localOpt.get());
    if (finalOpt) ransac.setFinalOptimizer(finalOpt.get());
    ransac.setSettings(settings);

    ransac.run(data);

    SuperansacStatistics stats;
    stats.num_iterations = ransac.getIterationNumber();
    stats.inliers = ransac.getInliers();

    PoseData best_pose;
    if (stats.inliers.size() >= estimator->sampleSize())
        best_pose = model_to_pose(ransac.getBestModel());

    return {std::move(best_pose), std::move(stats)};
}

} // namespace detail

template <class Solver>
std::pair<PoseData, SuperansacStatistics> run_superansac_point(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const superansac::RANSACSettings& settings)
{
    auto data = detail::pack_point_data(x, y);
    SuperansacPointEstimator<Solver> estimator(solver);
    return detail::run_pipeline<4>(&estimator, data, settings);
}

template <class Solver>
std::pair<PoseData, SuperansacStatistics> run_superansac_affine(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A,
    const superansac::RANSACSettings& settings)
{
    auto data = detail::pack_affine_data(x, y, A);
    SuperansacAffineEstimator<Solver> estimator(solver);
    return detail::run_pipeline<8>(&estimator, data, settings);
}

template <class Solver>
std::pair<PoseData, SuperansacStatistics> run_superansac_orientation(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori,
    const superansac::RANSACSettings& settings)
{
    auto data = detail::pack_orientation_data(x, y, ori);
    SuperansacOrientationEstimator<Solver> estimator(solver);
    return detail::run_pipeline<6>(&estimator, data, settings);
}

} // namespace HomLib

#endif // HOMLIB_HAS_SUPERANSAC
