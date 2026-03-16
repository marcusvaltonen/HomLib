#ifndef INCLUDES_HOMLIB_GET_WADENBACK_2025_HPP_
#define INCLUDES_HOMLIB_GET_WADENBACK_2025_HPP_

#include <Eigen/Dense>
#include <vector>
#include "posedata.hpp"
#include "radial.hpp"
#include "pose_estimator.h"
#include "refinement.hpp"
#include <iostream>

namespace HomLib {
namespace Wadenback2025 {
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
std::vector<HomLib::PoseData> get_unsided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
);
class SolverUnsided : public PoseEstimator<SolverUnsided> {
    public:
        SolverUnsided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback2025::get_unsided(x, y, extra_check);
            for (size_t i = 0; i < output.size(); i++) {
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 4;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            return xd;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            return yu;
        }        
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_unsided(x, y, pose);
        }
    private:
        bool extra_check = false;
    };
class SolverSingleSided : public PoseEstimator<SolverSingleSided> {
    public:
        SolverSingleSided() = default;
        int solve(const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y, std::vector<HomLib::PoseData> *poses) const {
            std::vector<HomLib::PoseData> output = HomLib::Wadenback2025::get_one_sided(x, y, extra_check);
            for (size_t i = 0; i < output.size(); i++) {
                // std::cout << "H[" << i << "] = " << output[i].homography / output[i].homography(2,2) << std::endl;
                // std::cout << "k1[" << i << "] = " << output[i].distortion_parameter << std::endl;
                // std::cout << "k2[" << i << "] = " << output[i].distortion_parameter2 << std::endl;
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, 0.0);  // One-sided
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter);
            return yd;
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
            std::vector<HomLib::PoseData> output = HomLib::Wadenback2025::get_double_sided_equal(x, y, extra_check);
            for (size_t i = 0; i < output.size(); i++) {
                // std::cout << "H[" << i << "] = " << output[i].homography / output[i].homography(2,2) << std::endl;
                // std::cout << "k1[" << i << "] = " << output[i].distortion_parameter << std::endl;
                // std::cout << "k2[" << i << "] = " << output[i].distortion_parameter2 << std::endl;
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter);
            return yd;
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
            std::vector<HomLib::PoseData> output = HomLib::Wadenback2025::get_double_sided(x, y, extra_check);
            for (size_t i = 0; i < output.size(); i++) {
                // std::cout << "H[" << i << "] = " << output[i].homography / output[i].homography(2,2) << std::endl;
                // std::cout << "k1[" << i << "] = " << output[i].distortion_parameter << std::endl;
                // std::cout << "k2[" << i << "] = " << output[i].distortion_parameter2 << std::endl;
                poses->push_back(output[i]);
            }
            return output.size();
        }
        int minimal_sample_size() const {
            return 5;
        }
        inline Eigen::Vector2d undistort(const HomLib::PoseData pose, const Eigen::Vector2d &xd) const {
            Eigen::Vector2d xu = HomLib::radialundistort(xd, pose.distortion_parameter);
            return xu;
        }
        inline Eigen::Vector2d distort(const HomLib::PoseData pose, const Eigen::Vector2d &yu) const {
            Eigen::Vector2d yd = HomLib::radialdistort(yu, pose.distortion_parameter2);
            return yd;
        }
        inline void refine(HomLib::PoseData &pose, const std::vector<Eigen::Vector2d> &x, const std::vector<Eigen::Vector2d> &y) const {
            HomLib::refinement_twosided(x, y, pose);
        }
    private:
        bool extra_check = false;
    };
}
}

#endif  // INCLUDES_HOMLIB_GET_WADENBACK_2025_HPP_
