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


#include <vector>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include "posedata.hpp"

#include "radial.hpp"
#include "get_nakano_icpr_2025.hpp"

namespace HomLib {
namespace NakanoICPR2025 {
    std::vector<HomLib::PoseData> get(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        bool extra_check
    ) {

        // Compute m vectors
        std::vector<Eigen::RowVector4d> m;
        for (int i = 0; i < 5; i++) {
            Eigen::RowVector4d tmp;
            tmp << y[i].squaredNorm(), y[i].transpose(), 1;
            m.push_back(tmp);
        }

        // Create M matrix
        Eigen::Matrix<double, 9, 12> M;
        M.setZero();
        for (int i = 0; i < 5; i++) {
            M.block(2*i,4,1,4) = -m[i];
            M.block(2*i,8,1,4) = x[i](1) * m[i];
            if (i < 4) {
                M.block(2*i+1,0,1,4) = m[i];
                M.block(2*i+1,8,1,4) = -x[i](0) * m[i];
            }
        }
        
        // Compute nullspace using QR
        Eigen::Matrix<double, 12, 12> Q = M.transpose().householderQr().householderQ();
        Eigen::Matrix<double, 12, 3> N = Q.rightCols(3);

        // Create generalized eigenvalue problem
        Eigen::Matrix3d A;
        A << N.row(0), N.row(4), N.row(8);
        Eigen::Matrix3d B;
        B << N.row(3), N.row(7), N.row(11);
        
        Eigen::GeneralizedEigenSolver<Eigen::Matrix3d> ges;
        ges.compute(A, B, true);
        Eigen::Vector3cd ks = ges.eigenvalues();
        Eigen::Matrix3d X = ges.eigenvectors().real();
        
        // Keep only real solutions (up to 3)
        std::vector<HomLib::PoseData> output;

        for (int i = 0; i < 3; i++) {
            if (std::abs(ks(i).imag()) < 1e-14) {
                double k = ks(i).real();
                Eigen::Vector3d alpha = X.col(i);
                Eigen::Matrix<double, 12, 1> g = alpha[0] * N.col(0) + alpha[1] * N.col(1) + alpha[2] * N.col(2);
                Eigen::Matrix<double, 3, 4> G = g.reshaped(4, 3).transpose();
                Eigen::Matrix3d H = G.bottomRightCorner(3, 3);
                
                // Package output
                HomLib::PoseData pd;
                pd.homography = H.inverse();
                pd.distortion_parameter = k;
                output.push_back(pd);
            }
        }
        
        // Compute reprojection error using the unused constraint
        if (extra_check) {
            double min_res = std::numeric_limits<double>::max();
            int best_id = -1;
            for (size_t i = 0; i < output.size(); i++) {
                // Measure reprojection error in undistorted space
                Eigen::Vector2d y4u_est1 = HomLib::radialundistort(y[4], output[i].distortion_parameter);
                Eigen::Vector3d tmp = output[i].homography * x[4].homogeneous();
                Eigen::Vector2d y4u_est2 = tmp.hnormalized();
                double res = (y4u_est1 - y4u_est2).squaredNorm();
                if (res < min_res) {
                    min_res = res;
                    best_id = i;
                }
            }
            std::vector<HomLib::PoseData> output2;
            if (best_id >= 0) {
                output2.push_back(output[best_id]);
            }
            return output2;
        }
        return output;
    }
}  // namespace NakanoICPR2025
}  // namespace HomLib
