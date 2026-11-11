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
#include "neighborhood/grid_neighborhood_graph.h"
#include "neighborhood/flann_neighborhood_graph.h"
#include "local_optimization/graph_cut_ransac_optimizer.h"
#include "samplers/napsac_sampler.h"
#include "samplers/progressive_napsac_sampler.h"
#include "samplers/importance_sampler.h"
#include "samplers/adaptive_reordering_sampler.h"
#include "scoring/grid_scoring.h"
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

template <size_t DimensionNumber>
void initialize_neighborhood(
    const DataMatrix& data,
    std::unique_ptr<superansac::neighborhood::AbstractNeighborhoodGraph>& graph,
    const superansac::RANSACSettings& settings)
{
    using namespace superansac::neighborhood;
    graph = createNeighborhoodGraph<DimensionNumber>(settings.neighborhood);

    if (settings.neighborhood == NeighborhoodType::Grid) {
        auto* grid = dynamic_cast<GridNeighborhoodGraph<DimensionNumber>*>(graph.get());
        const double density = settings.neighborhoodSettings.neighborhoodGridDensity;
        std::vector<double> cellSizes(DimensionNumber);
        for (size_t i = 0; i < DimensionNumber; ++i) {
            double range = data.col(i).maxCoeff() - data.col(i).minCoeff();
            cellSizes[i] = std::max(1.0, range / density);
        }
        grid->setCellSizes(cellSizes, static_cast<size_t>(density));
    } else if (settings.neighborhood == NeighborhoodType::FLANN_KNN) {
        auto* flann = dynamic_cast<FlannNeighborhoodGraph<DimensionNumber, 0>*>(graph.get());
        flann->setNearestNeighborNumber(
            settings.neighborhoodSettings.nearestNeighborNumber);
    } else if (settings.neighborhood == NeighborhoodType::FLANN_Radius) {
        auto* flann = dynamic_cast<FlannNeighborhoodGraph<DimensionNumber, 1>*>(graph.get());
        flann->setRadius(settings.neighborhoodSettings.neighborhoodSize);
    }
    graph->initialize(&data);
}

inline bool needs_neighborhood(
    superansac::local_optimization::LocalOptimizationType lo,
    superansac::local_optimization::LocalOptimizationType fo,
    superansac::samplers::SamplerType sampler)
{
    using LO = superansac::local_optimization::LocalOptimizationType;
    using S = superansac::samplers::SamplerType;
    return lo == LO::GCRANSAC || fo == LO::GCRANSAC ||
           sampler == S::NAPSAC || sampler == S::ProgressiveNAPSAC;
}

template <size_t DimensionNumber>
void initialize_local_optimizer(
    std::unique_ptr<superansac::local_optimization::LocalOptimizer>& opt,
    superansac::local_optimization::LocalOptimizationType type,
    const superansac::LocalOptimizationSettings& loSettings,
    bool isFinal,
    superansac::neighborhood::AbstractNeighborhoodGraph* neighborhood = nullptr)
{
    using namespace superansac::local_optimization;
    switch (type) {
        case LocalOptimizationType::IRLS: {
            auto* irls = dynamic_cast<IRLSOptimizer*>(opt.get());
            irls->setMaxIterations(loSettings.maxIterations);
            irls->setUseInliers(true);
            break;
        }
        case LocalOptimizationType::LSQ: {
            auto* lsq = dynamic_cast<LeastSquaresOptimizer*>(opt.get());
            lsq->setUseInliers(true);
            break;
        }
        case LocalOptimizationType::NestedRANSAC: {
            auto* nested = dynamic_cast<NestedRANSACOptimizer*>(opt.get());
            nested->setMaxIterations(loSettings.maxIterations);
            nested->setSampleSizeMultiplier(loSettings.sampleSizeMultiplier);
            break;
        }
        case LocalOptimizationType::GCRANSAC: {
            auto* gc = dynamic_cast<GraphCutRANSACOptimizer*>(opt.get());
            if (neighborhood)
                gc->setNeighborhood(neighborhood);
            gc->setMaxIterations(loSettings.maxIterations);
            gc->setGraphCutNumber(loSettings.graphCutNumber);
            gc->setSampleSizeMultiplier(loSettings.sampleSizeMultiplier);
            gc->setSpatialCoherenceWeight(loSettings.spatialCoherenceWeight);
            break;
        }
        default:
            break;
    }
}

// SupeRANSAC's UniformRandomSampler, initialized on the data, draws indices in [0, rows], one past
// the last point (danini/superansac#10): this one draws them in [0, rows - 1]. It sets the range
// itself, so it stays right whichever fix SupeRANSAC adopts (#9 also changes initialize(size_t),
// which its local optimizers call with the largest index).
class UniformSampler : public superansac::samplers::UniformRandomSampler {
public:
    void initialize(const DataMatrix& data) override {
        randomGenerator = std::make_unique<superansac::utils::UniformRandomGenerator<size_t>>();
        randomGenerator->resetGenerator(0, data.rows() > 0 ? static_cast<size_t>(data.rows()) - 1 : 0);
    }
};

// The sampler of the given type, with HomLib's uniform sampler instead of SupeRANSAC's
template <size_t DimensionNumber>
std::unique_ptr<superansac::samplers::AbstractSampler> create_sampler(superansac::samplers::SamplerType type)
{
    if (type == superansac::samplers::SamplerType::Uniform)
        return std::make_unique<UniformSampler>();
    return superansac::samplers::createSampler<DimensionNumber>(type);
}

// Assemble and run the SupeRANSAC pipeline with a given estimator and data.
// DimensionNumber is the number of columns in the DataMatrix (4, 6, or 8).
template <size_t DimensionNumber>
std::pair<PoseData, SuperansacStatistics> run_pipeline(
    superansac::estimator::Estimator* estimator,
    const DataMatrix& data,
    const std::vector<double>& probabilities,
    const superansac::RANSACSettings& settings)
{
    // Neighborhood graph (lazily initialized for GCRANSAC LO, NAPSAC samplers, GRID scoring)
    std::unique_ptr<superansac::neighborhood::AbstractNeighborhoodGraph> neighborhoodGraph;
    bool useNeighborhood = needs_neighborhood(
        settings.localOptimization, settings.finalOptimization, settings.sampler) ||
        settings.scoring == superansac::scoring::ScoringType::GRID;

    if (useNeighborhood)
        initialize_neighborhood<DimensionNumber>(data, neighborhoodGraph, settings);

    auto sampler = create_sampler<DimensionNumber>(settings.sampler);
    auto scoring = superansac::scoring::createScoring<DimensionNumber>(settings.scoring, settings.useSprt);
    auto termination = superansac::termination::createTerminationCriterion(settings.terminationCriterion);

    // Sampler initialization
    using S = superansac::samplers::SamplerType;
    if (settings.sampler == S::PROSAC)
        dynamic_cast<superansac::samplers::PROSACSampler*>(sampler.get())
            ->setSampleSize(estimator->sampleSize());
    if (settings.sampler == S::NAPSAC && neighborhoodGraph)
        dynamic_cast<superansac::samplers::NAPSACSampler*>(sampler.get())
            ->setNeighborhood(neighborhoodGraph.get());
    if (!probabilities.empty() && settings.sampler == S::ARSampler) {
        dynamic_cast<superansac::samplers::AdaptiveReorderingSampler*>(sampler.get())
            ->initialize(&data, probabilities);
    } else if (!probabilities.empty() && settings.sampler == S::ImportanceSampler) {
        dynamic_cast<superansac::samplers::ImportanceSampler*>(sampler.get())
            ->setProbabilities(probabilities);
        sampler->initialize(data);
    } else {
        sampler->initialize(data);
    }

    // Scoring initialization
    scoring->setThreshold(settings.inlierThreshold);
    if (settings.scoring == superansac::scoring::ScoringType::MAGSAC) {
        if (settings.useSprt)
            dynamic_cast<superansac::scoring::MAGSACSPRTScoring*>(scoring.get())->initialize(estimator);
        else
            dynamic_cast<superansac::scoring::MAGSACScoring*>(scoring.get())->initialize(estimator);
    }
    if (settings.scoring == superansac::scoring::ScoringType::GRID && neighborhoodGraph) {
        auto* gridGraph = dynamic_cast<superansac::neighborhood::GridNeighborhoodGraph<DimensionNumber>*>(
            neighborhoodGraph.get());
        if (gridGraph)
            dynamic_cast<superansac::scoring::GridScoring<DimensionNumber>*>(scoring.get())
                ->setNeighborhood(gridGraph);
    }

    // Termination criterion initialization
    if (settings.terminationCriterion == superansac::termination::TerminationType::RANSAC)
        dynamic_cast<superansac::termination::RANSACCriterion*>(termination.get())
            ->setConfidence(settings.confidence);

    // Local and final optimizer initialization
    superansac::neighborhood::AbstractNeighborhoodGraph* nbPtr =
        neighborhoodGraph ? neighborhoodGraph.get() : nullptr;

    std::unique_ptr<superansac::local_optimization::LocalOptimizer> localOpt, finalOpt;
    if (settings.localOptimization != superansac::local_optimization::LocalOptimizationType::None) {
        localOpt = superansac::local_optimization::createLocalOptimizer(settings.localOptimization);
        initialize_local_optimizer<DimensionNumber>(
            localOpt, settings.localOptimization, settings.localOptimizationSettings, false, nbPtr);
    }
    if (settings.finalOptimization != superansac::local_optimization::LocalOptimizationType::None) {
        finalOpt = superansac::local_optimization::createLocalOptimizer(settings.finalOptimization);
        initialize_local_optimizer<DimensionNumber>(
            finalOpt, settings.finalOptimization, settings.finalOptimizationSettings, true, nbPtr);
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
    const std::vector<double>& probabilities,
    const NonMinimalSettings& nm,
    const superansac::RANSACSettings& settings)
{
    auto data = detail::pack_point_data(x, y);
    SuperansacPointEstimator<Solver> estimator(solver, nm);
    auto result = detail::run_pipeline<4>(&estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

template <class Solver>
std::pair<PoseData, SuperansacStatistics> run_superansac_affine(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A,
    const std::vector<double>& probabilities,
    const NonMinimalSettings& nm,
    const superansac::RANSACSettings& settings)
{
    auto data = detail::pack_affine_data(x, y, A);
    SuperansacAffineEstimator<Solver> estimator(solver, nm);
    auto result = detail::run_pipeline<8>(&estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

template <class Solver>
std::pair<PoseData, SuperansacStatistics> run_superansac_orientation(
    const Solver& solver,
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori,
    const std::vector<double>& probabilities,
    const NonMinimalSettings& nm,
    const superansac::RANSACSettings& settings)
{
    auto data = detail::pack_orientation_data(x, y, ori);
    SuperansacOrientationEstimator<Solver> estimator(solver, nm);
    auto result = detail::run_pipeline<6>(&estimator, data, probabilities, settings);
    detail::refine_on_inliers(solver, x, y, result.second.inliers, result.first);
    return result;
}

} // namespace HomLib

#endif // HOMLIB_HAS_SUPERANSAC
