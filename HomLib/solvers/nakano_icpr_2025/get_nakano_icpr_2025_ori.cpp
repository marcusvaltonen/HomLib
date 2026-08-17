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
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include "posedata.hpp"

#include "radial.hpp"
#include "get_nakano_icpr_2025.hpp"

namespace HomLib {
namespace NakanoICPR2025 {
    std::vector<HomLib::PoseData> get_ori(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Vector2d> &ori
    ) {

        double weight = 1.0;

        // Create M matrix
        Eigen::Matrix<double, Eigen::Dynamic, 12> M(3 * x.size(), 12);

        size_t rowIdx = 0;
        for (size_t i = 0; i < x.size(); i++) {
            double
	            x1 = x[i](0),
	            y1 = x[i](1),
	            x2 = y[i](0),
	            y2 = y[i](1),
	            s1 = std::sin(ori[i](0)),
	            s2 = std::sin(ori[i](1)),
	            c1 = std::cos(ori[i](0)),
	            c2 = std::cos(ori[i](1));

	        double r12 = x1 * x1 + y1 * y1;

            double
                kMinusWeightTimesR12 = -weight * r12,
	            kMinusWeightTimesX1 = -weight * x1,
	            kMinusWeightTimesY1 = -weight * y1,
	            kWeightTimesX2 = weight * x2,
	            kWeightTimesY2 = weight * y2;

            M(rowIdx, 0) = kMinusWeightTimesR12;
            M(rowIdx, 1) = kMinusWeightTimesX1;
            M(rowIdx, 2) = kMinusWeightTimesY1;
            M(rowIdx, 3) = -weight;
            M(rowIdx, 4) = 0;
            M(rowIdx, 5) = 0;
            M(rowIdx, 6) = 0;
            M(rowIdx, 7) = 0;
            M(rowIdx, 8) = kWeightTimesX2 * r12;
            M(rowIdx, 9) = kWeightTimesX2 * x1;
            M(rowIdx, 10) = kWeightTimesX2 * y1;
            M(rowIdx, 11) = kWeightTimesX2;
            ++rowIdx;

            M(rowIdx, 0) = 0;
            M(rowIdx, 1) = 0;
            M(rowIdx, 2) = 0;
            M(rowIdx, 3) = 0;
            M(rowIdx, 4) = kMinusWeightTimesR12;
            M(rowIdx, 5) = kMinusWeightTimesX1;
            M(rowIdx, 6) = kMinusWeightTimesY1;
            M(rowIdx, 7) = -weight;
            M(rowIdx, 8) = kWeightTimesY2 * r12;
            M(rowIdx, 9) = kWeightTimesY2 * x1;
            M(rowIdx, 10) = kWeightTimesY2 * y1;
            M(rowIdx, 11) = kWeightTimesY2;
            ++rowIdx;

            M(rowIdx, 0) = -2 * (x1*s2*c1 + y1*s1*s2);
            M(rowIdx, 1) = -s2 * c1;
            M(rowIdx, 2) = -s1 * s2;
            M(rowIdx, 3) = 0;
            M(rowIdx, 4) = 2 * (x1*c1*c2 + y1*s1*c2);
            M(rowIdx, 5) = c1*c2;
            M(rowIdx, 6) = s1*c2;
            M(rowIdx, 7) = 0;
            M(rowIdx, 8) = 2*x1*(x2*s2*c1 - y2*c1*c2) + 2*y1*(x2*s1*s2 - y2*s1*c2);
            M(rowIdx, 9) = x2*s2*c1 - y2*c1*c2;
            M(rowIdx, 10) = x2*s1*s2 - y2*s1*c2;
            M(rowIdx, 11) = 0;
            ++rowIdx;
        }

        // Compute nullspace using QR
        Eigen::Matrix<double, 12, 12> Q = M.transpose().householderQr().householderQ();
        
        std::vector<HomLib::PoseData> output;
        if (M.rows() < 12) {
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
        } else {
            // Non-minimal
            Eigen::Matrix<double, 12, 1> g = Q.rightCols(1);
            Eigen::Matrix<double, 3, 4> G = g.reshaped(4, 3).transpose();
            Eigen::Matrix3d H = G.bottomRightCorner(3, 3);
            // Package output
            HomLib::PoseData pd;
            pd.homography = H;  // Note: not inverse now
            pd.distortion_parameter = (double) (G.col(0).transpose() * G.col(3)) / G.col(3).squaredNorm();
            pd.distortion_parameter2 = 0.0;
            pd.focal_length = 0.0;
            output.push_back(pd);
        }
        return output;
    }

}  // namespace NakanoICPR2025
}  // namespace HomLib
