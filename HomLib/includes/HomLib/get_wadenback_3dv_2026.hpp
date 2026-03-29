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
std::vector<HomLib::PoseData> get_double_sided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
);
std::vector<HomLib::PoseData> get_double_sided_equal(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
);
std::vector<HomLib::PoseData> get_one_sided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
);
class SolverSingleSided : public PoseEstimator<SolverSingleSided> {
    public:
        SolverSingleSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback3DV2026::get_one_sided(x, y, extra_check);
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
        SolverTwoSidedEqual() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback3DV2026::get_double_sided_equal(x, y, extra_check);
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
        SolverTwoSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback3DV2026::get_double_sided(x, y, extra_check);
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
