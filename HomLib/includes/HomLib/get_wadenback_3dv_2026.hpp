#ifndef INCLUDES_HOMLIB_GET_WADENBACK_3DV_2026_HPP_
#define INCLUDES_HOMLIB_GET_WADENBACK_3DV_2026_HPP_

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "pose_estimator.h"
#include "refinement.hpp"
#include <iostream>

namespace HomLib {
namespace Wadenback3DV2026 {
// With a range: keeping only the solutions with distortion parameters in the range (checked before
// computing their homography, and before choosing the best root with extra_check)
std::vector<HomLib::PoseData> get_double_sided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check,
    const HomLib::DistortionRange &range = HomLib::DistortionRange()
);
// With a range: keeping only the solutions with distortion parameters in the range (checked before
// computing their homography, and before choosing the best root with extra_check)
std::vector<HomLib::PoseData> get_double_sided_equal(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check,
    const HomLib::DistortionRange &range = HomLib::DistortionRange()
);
// With a range: keeping only the solutions with distortion parameters in the range (checked before
// computing their homography, and before choosing the best root with extra_check)
std::vector<HomLib::PoseData> get_one_sided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check,
    const HomLib::DistortionRange &range = HomLib::DistortionRange()
);
class SolverSingleSided : public PoseEstimator<SolverSingleSided> {
    public:
        static constexpr int minimal_sample_size_value = 5;
        SolverSingleSided() = default;
        explicit SolverSingleSided(bool extra_check) : extra_check(extra_check) {}
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback3DV2026::get_one_sided(x, y, extra_check, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_onesided(x, y, pose);
        }
    private:
        bool extra_check = true;
    };
class SolverTwoSidedEqual : public PoseEstimator<SolverTwoSidedEqual> {
    public:
        static constexpr int minimal_sample_size_value = 5;
        SolverTwoSidedEqual() = default;
        explicit SolverTwoSidedEqual(bool extra_check) : extra_check(extra_check) {}
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback3DV2026::get_double_sided_equal(x, y, extra_check, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_twosided_equal(x, y, pose);
        }
    private:
        bool extra_check = true;
    };
class SolverTwoSided : public PoseEstimator<SolverTwoSided> {
    public:
        static constexpr int minimal_sample_size_value = 5;
        // With the range check of the estimator (PoseEstimator::use_k_range) on by default: of the
        // 9 roots of the solver (about 7.8 real ones per noise-free normalized sample), most have
        // implausible distortion parameters (about 1.8 are kept), and scoring fewer hypotheses makes
        // LO-MSAC 1.6 to 2.2 times faster at the same accuracy. This check isn't part of the paper's
        // method. The default range, [-1, 0.1] in the normalized coordinates, covers real cameras
        // (about [-0.45, 0.04] there, derived from the distortion measured on hundreds of cameras by
        // Lopez-Antequera et al., CVPR 2019, Fig. 3) with a margin for noise.
        SolverTwoSided() { use_k_range = true; }
        explicit SolverTwoSided(bool extra_check) : extra_check(extra_check) { use_k_range = true; }
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses,
                  const HomLib::DistortionRange &range = HomLib::DistortionRange()) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback3DV2026::get_double_sided(x, y, extra_check, range);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_twosided(x, y, pose);
        }
    private:
        bool extra_check = false;
    };
}
}

#endif  // INCLUDES_HOMLIB_GET_WADENBACK_3DV_2026_HPP_
