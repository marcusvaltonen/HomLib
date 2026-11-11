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

#include "invert_homography.hpp"
#include "radial.hpp"
#include "get_nakano_icpr_2025.hpp"

namespace HomLib {
namespace NakanoICPR2025 {
    std::vector<HomLib::PoseData> get(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        bool extra_check,
        const HomLib::DistortionRange &range
    ) {
        const int n_points = static_cast<int>(x.size());
        if (n_points < 5)
            return {};

        // Compute m vectors
        std::vector<Eigen::RowVector4d> m;
        for (int i = 0; i < n_points; i++) {
            Eigen::RowVector4d tmp;
            tmp << y[i].squaredNorm(), y[i].transpose(), 1;
            m.push_back(tmp);
        }

        // Compute 3-dim nullspace N.
        // Minimal (n_points == 5): fixed-size 9x12 matrix, QR.
        // Overdetermined (n_points > 5): SVD of M (not M^T M, which squares its condition).
        Eigen::Matrix<double, 12, 3> N;
        if (n_points == 5) {
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
            Eigen::Matrix<double, 12, 12> Q = M.transpose().householderQr().householderQ();
            N = Q.rightCols(3);
        } else {
            Eigen::MatrixXd M(2 * n_points, 12);
            M.setZero();
            for (int i = 0; i < n_points; i++) {
                M.block(2*i,4,1,4) = -m[i];
                M.block(2*i,8,1,4) = x[i](1) * m[i];
                M.block(2*i+1,0,1,4) = m[i];
                M.block(2*i+1,8,1,4) = -x[i](0) * m[i];
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullV);
            N = svd.matrixV().rightCols(3);
        }

        // Create generalized eigenvalue problem
        Eigen::Matrix3d A;
        A << N.row(0), N.row(4), N.row(8);
        Eigen::Matrix3d B;
        B << N.row(3), N.row(7), N.row(11);

        // Eigen's RealQZ doesn't terminate when A = 0 (and B is singular), as with all the points at
        // the origin
        if (A.isZero(0.0)) {
            return {};
        }

        Eigen::GeneralizedEigenSolver<Eigen::Matrix3d> ges;
        ges.compute(A, B, true);
        if (ges.info() != Eigen::Success) {
            return {};
        }
        Eigen::Vector3cd ks = ges.eigenvalues();
        Eigen::Matrix3d X = ges.eigenvectors().real();

        // Keep only real solutions (up to 3)
        std::vector<HomLib::PoseData> output;

        for (int i = 0; i < 3; i++) {
            if (std::abs(ks(i).imag()) < 1e-14) {
                double k = ks(i).real();
                if (!range.contains(k)) {
                    continue;
                }
                Eigen::Vector3d alpha = X.col(i);
                Eigen::Matrix<double, 12, 1> g = alpha[0] * N.col(0) + alpha[1] * N.col(1) + alpha[2] * N.col(2);
                Eigen::Matrix<double, 3, 4> G = Eigen::Map<const Eigen::Matrix<double, 4, 3>>(g.data()).transpose();
                Eigen::Matrix3d H = G.bottomRightCorner(3, 3);

                // Package output; a singular H gives no solution
                HomLib::PoseData pd;
                if (!HomLib::invert_homography(H, &pd.homography))
                    continue;
                pd.distortion_parameter2 = k;
                output.push_back(pd);
            }
        }

        // Select best solution using reprojection error
        if (extra_check && !output.empty()) {
            double min_res = std::numeric_limits<double>::max();
            int best_id = -1;
            for (size_t i = 0; i < output.size(); i++) {
                double res = 0.0;
                for (int j = 0; j < n_points; j++) {
                    Eigen::Vector2d yu_est1;
                    if (!HomLib::try_radialundistort(y[j], output[i].distortion_parameter2, yu_est1)) {
                        res = std::numeric_limits<double>::max();
                        break;
                    }
                    Eigen::Vector3d tmp = output[i].homography * x[j].homogeneous();
                    Eigen::Vector2d yu_est2 = tmp.hnormalized();
                    res += (yu_est1 - yu_est2).squaredNorm();
                }
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
