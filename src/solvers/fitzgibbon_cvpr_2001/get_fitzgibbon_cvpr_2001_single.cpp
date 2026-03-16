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

#include "get_fitzgibbon_cvpr_2001.hpp"
#include <vector>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include "posedata.hpp"
#include "radial.hpp"

namespace HomLib {
namespace FitzgibbonCVPR2001 {
    std::vector<HomLib::PoseData> get_single(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y
    ) {
        // This is a five point method
        const int n_points = 5;

        // Compute the distance to center point
        std::vector<double> r2;
        for (int i = 0; i < n_points; i++) {
            r2.push_back(y[i].squaredNorm());
        }

        // Follow Nakano's nomenclature
        Eigen::Matrix<double, 9, 9> D0;
        D0.setZero();
        Eigen::Matrix<double, 9, 9> D1;
        D1.setZero();

        for (int i = 0; i < n_points; i++) {
            D0.block(2*i,3,1,3) = -x[i].homogeneous().transpose();
            D0.block(2*i,6,1,3) = y[i](1) * x[i].homogeneous().transpose();
            if (i < 4) {
                D0.block(2*i+1,0,1,3) = x[i].homogeneous().transpose();
                D0.block(2*i+1,6,1,3) = -y[i](0) * x[i].homogeneous().transpose();
            }
            
            D1.block(2*i,3,1,3) = r2[i]*x[i].homogeneous().transpose();
            if (i < 4) {
                D1.block(2*i+1,0,1,3) = -r2[i]*x[i].homogeneous().transpose();
            }
        }

        // Create and solve generalized eigenvalue problem
        // In the minimal case, we don't need to compute D0.transpose * D0, etc.
        Eigen::GeneralizedEigenSolver<Eigen::Matrix<double, 9, 9>> ges;
        ges.compute(D0, D1, true);
        Eigen::Matrix<std::complex<double>, 9, 1> ks = ges.eigenvalues();
        Eigen::Matrix<double, 9, 9> X = ges.eigenvectors().real();

        // Extract correct solution
        double min_res = std::numeric_limits<double>::max();
        std::vector<HomLib::PoseData> output;
        HomLib::PoseData posedata;

        for (int i = 0; i < 9; i++) {
            if (std::abs(ks(i).imag()) < 1e-14) {
                double k = ks(i).real();
                Eigen::Matrix3d H = Eigen::Map<Eigen::Matrix3d>(X.col(i).data(), 3, 3).transpose();
                // Measure reprojection error in undistorted space
                Eigen::Vector2d y4u_est1 = HomLib::radialundistort(y[4], k);
                Eigen::Vector3d tmp = H * x[4].homogeneous();
                Eigen::Vector2d y4u_est2 = tmp.hnormalized();
                double res = (y4u_est1 - y4u_est2).squaredNorm();
                if (res < min_res) {
                    min_res = res;
                    posedata.homography = H;
                    posedata.distortion_parameter = k;
                }
            }
        }
        if (min_res < std::numeric_limits<double>::max())
            output.push_back(posedata);

        return output;
    }
}  // namespace FitzgibbonCVPR2001
}  // namespace HomLib
