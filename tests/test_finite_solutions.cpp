#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "generate_problem_instance.hpp"
#include "posedata.hpp"
#include "problem_instance.hpp"
#include "solver_dispatch.h"

// Every solver, through estimate() (normalized, with the robust estimators' range) and solve()
// (as the Python bindings call it), on synthetic samples and degenerate ones, returns only finite
// solutions (homographies and distortion parameters), and its refinement keeps them finite.

namespace {

using V2 = Eigen::Vector2d;
using Poses = std::vector<HomLib::PoseData>;
using Affines = std::vector<Eigen::Matrix2d>;

bool finite(const HomLib::PoseData &p) {
    return p.homography.allFinite() && std::isfinite(p.distortion_parameter) &&
           std::isfinite(p.distortion_parameter2) && std::isfinite(p.focal_length);
}

// estimate() or solve(), with the arguments of the solver's feature type
template <class S>
void call(const S &s, bool direct, const std::vector<V2> &x, const std::vector<V2> &y, const Affines &A,
          const std::vector<V2> &ori, Poses *poses) {
    const HomLib::DistortionRange none;
    if constexpr (std::is_base_of_v<HomLib::AffinePoseEstimator<S>, S>)
        direct ? (void)s.solve(x, y, A, poses, none) : (void)s.estimate(x, y, A, poses);
    else if constexpr (std::is_base_of_v<HomLib::OrientationPoseEstimator<S>, S>)
        direct ? (void)s.solve(x, y, ori, poses, none) : (void)s.estimate(x, y, ori, poses);
    else
        direct ? (void)s.solve(x, y, poses, none) : (void)s.estimate(x, y, poses);
}

} // anonymous namespace

TEST_CASE("Solvers - every returned solution is finite, also on degenerate samples") {
    // Every SolverId: dispatch_solver throws for the first value past the last one
    int n_solvers = 0;
    for (int id = 0;; ++id, ++n_solvers) {
        try {
            HomLib::dispatch_solver(HomLib::SolverId(id), [&](auto solver) {
                const int n = decltype(solver)::minimal_sample_size_value;
                for (int dc = HomLib::NO_DISTORTION; dc <= HomLib::TWO_SIDED; ++dc)
                    for (int seed = 0; seed < 10; ++seed) {
                        HomLib::ProblemConfig config;
                        config.distortion = HomLib::DistortionCase(dc);
                        config.point_noise = seed % 2 ? 1e-3 : 0.0;
                        config.number_points = 4 * n;
                        config.seed = seed;
                        const HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);
                        const std::vector<V2> x(inst.x1.begin(), inst.x1.begin() + n), y(inst.x2.begin(), inst.x2.begin() + n);
                        const Affines A(inst.A.begin(), inst.A.begin() + n);
                        const std::vector<V2> ori(inst.ori.begin(), inst.ori.begin() + n);

                        std::vector<std::pair<std::string, std::pair<std::vector<V2>, std::vector<V2>>>> samples;
                        samples.push_back({"synthetic", {x, y}});
                        {
                            auto x2 = x, y2 = y;
                            x2[n - 1] = x2[0];
                            y2[n - 1] = y2[0];
                            samples.push_back({"repeated point", {x2, y2}});
                        }
                        {
                            auto x2 = x, y2 = y;
                            for (int i = 0; i < n; ++i) {
                                x2[i] = x[0] + (i - 1.5) * (x[1] - x[0]);
                                y2[i] = y[0] + (i - 1.5) * (y[1] - y[0]);
                            }
                            samples.push_back({"collinear", {x2, y2}});
                        }
                        {
                            auto x2 = x, y2 = y;
                            x2[0].setZero();
                            y2[0].setZero();
                            samples.push_back({"a point at the origin", {x2, y2}});
                        }
                        samples.push_back({"all at the origin", {std::vector<V2>(n, V2::Zero()), std::vector<V2>(n, V2::Zero())}});
                        samples.push_back({"y = x", {x, x}});
                        {
                            auto x2 = x, y2 = y;
                            for (int i = 0; i < n; ++i) {
                                x2[i] *= 500;   // pixels
                                y2[i] *= 500;
                            }
                            samples.push_back({"scaled by 500", {x2, y2}});
                        }

                        for (const auto &[what, xy] : samples)
                            for (bool direct : {false, true}) {
                                INFO("solver " << id << ", " << what << (direct ? ", solve()" : ", estimate()")
                                     << ", distortion case " << dc << ", seed " << seed);
                                Poses poses;
                                call(solver, direct, xy.first, xy.second, A, ori, &poses);
                                for (auto pose : poses) {
                                    CHECK(finite(pose));
                                    if (!finite(pose))
                                        continue;
                                    solver.refine(pose, inst.x1, inst.x2);
                                    CHECK(finite(pose));
                                }
                            }
                    }
            });
        } catch (const std::invalid_argument &) {
            break;
        }
    }
    CHECK(n_solvers >= 17);
}
