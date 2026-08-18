// Copyright (c) 2020 Marcus Valtonen Örnhag
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef SRC_HELPERS_PROBLEM_INSTANCE_HPP_
#define SRC_HELPERS_PROBLEM_INSTANCE_HPP_

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include "posedata.hpp"

namespace HomLib {

enum DistortionCase {
    NO_DISTORTION,
    ONE_SIDED_LEFT,
    ONE_SIDED_RIGHT,
    TWO_SIDED_EQUAL,
    TWO_SIDED
}; 

struct ProblemInstance {
    HomLib::PoseData posedata;
    std::vector<Eigen::Vector2d> x1;
    std::vector<Eigen::Vector2d> x2;
    std::vector<Eigen::Matrix2d> A;
    std::vector<Eigen::Vector2d> ori;

    double hom_error(const Eigen::Matrix3d &H_est) const
    {   
        Eigen::Matrix3d H1 = posedata.homography.normalized();
        if (H1(2,2) < 0)
            H1 *= -1.0;
        Eigen::Matrix3d H2 = H_est.normalized();
        if (H2(2,2) < 0)
            H2 *= -1.0;
        return (H1 - H2).norm();
    }
    double dist_error(double k1_est, double k2_est) const
    {   
        return 0.5*(std::abs(posedata.distortion_parameter-k1_est) + std::abs(posedata.distortion_parameter2-k2_est));
    }
};

struct ProblemConfig {
    DistortionCase distortion;
    double point_noise;
    int number_points;
    double camera_fov_ = 70.0;
    double min_depth_ = 0.1;
    double max_depth_ = 10.0;
    double min_focal_ = 100.0;
    double max_focal_ = 1000.0;
    double min_dist_ = -0.2;
    double max_dist_ = -0.01;
    bool no_distortion = false;
};
}

#endif  // SRC_HELPERS_PROBLEM_INSTANCE_HPP_
