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
        const std::vector<Eigen::Vector2d> &y,
        const HomLib::DistortionRange &range
    ) {
        const int n_points = static_cast<int>(x.size());
        if (n_points < 5)
            return {};

        std::vector<double> r2(n_points);
        for (int i = 0; i < n_points; i++) {
            r2[i] = y[i].squaredNorm();
        }

        // Build 9x9 GEP matrices A, B.
        // Minimal (n_points == 5): fixed-size 9x9, A = D0, B = D1.
        // Overdetermined (n_points > 5): QR-project the pencil onto D0 (avoids squaring its
        // condition number). Still biased with noise; refine the result.
        Eigen::Matrix<double, 9, 9> A, B;
        if (n_points == 5) {
            Eigen::Matrix<double, 9, 9> D0, D1;
            D0.setZero();
            D1.setZero();
            for (int i = 0; i < 5; i++) {
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
            A = D0;
            B = D1;
        } else {
            const int n_rows = 2 * n_points;
            Eigen::MatrixXd D0(n_rows, 9), D1(n_rows, 9);
            D0.setZero();
            D1.setZero();
            for (int i = 0; i < n_points; i++) {
                D0.block(2*i,3,1,3) = -x[i].homogeneous().transpose();
                D0.block(2*i,6,1,3) = y[i](1) * x[i].homogeneous().transpose();
                D0.block(2*i+1,0,1,3) = x[i].homogeneous().transpose();
                D0.block(2*i+1,6,1,3) = -y[i](0) * x[i].homogeneous().transpose();
                D1.block(2*i,3,1,3) = r2[i]*x[i].homogeneous().transpose();
                D1.block(2*i+1,0,1,3) = -r2[i]*x[i].homogeneous().transpose();
            }
            // QR projection: avoids squaring D0's condition number (normal
            // equations fail when k ≈ 0 because D0 becomes rank-deficient).
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(D0);
            Eigen::MatrixXd thin_Q = qr.householderQ()
                * Eigen::MatrixXd::Identity(n_rows, 9);
            A.setZero();
            A.triangularView<Eigen::Upper>() = qr.matrixQR().topRows(9);
            B = thin_Q.transpose() * D1;
        }

        Eigen::GeneralizedEigenSolver<Eigen::Matrix<double, 9, 9>> ges;
        ges.compute(A, B, true);
        if (ges.info() != Eigen::Success) {
            return {};
        }
        Eigen::Matrix<std::complex<double>, 9, 1> ks = ges.eigenvalues();
        Eigen::Matrix<double, 9, 9> X = ges.eigenvectors().real();

        double min_res = std::numeric_limits<double>::max();
        std::vector<HomLib::PoseData> output;
        HomLib::PoseData posedata;

        for (int i = 0; i < 9; i++) {
            if (std::abs(ks(i).imag()) < 1e-14) {
                double k = ks(i).real();
                if (!range.contains(k)) {
                    continue;
                }
                Eigen::Matrix3d H = Eigen::Map<Eigen::Matrix3d>(X.col(i).data(), 3, 3).transpose();
                double res = 0.0;
                for (int j = 0; j < n_points; j++) {
                    Eigen::Vector2d yu_j;
                    if (!HomLib::try_radialundistort(y[j], k, yu_j)) {
                        res = std::numeric_limits<double>::max();
                        break;
                    }
                    Eigen::Vector3d tmp = H * x[j].homogeneous();
                    Eigen::Vector2d yu_est_j = tmp.hnormalized();
                    res += (yu_j - yu_est_j).squaredNorm();
                }
                if (res < min_res) {
                    min_res = res;
                    posedata.homography = H;
                    posedata.distortion_parameter2 = k;
                }
            }
        }
        // k = 0 fallback: the GEP is singular there, so recover it from the
        // nullspace of D0 (the standard DLT solution, no distortion).
        if (range.contains(0.0)) {
            Eigen::MatrixXd D0_full(2 * n_points, 9);
            D0_full.setZero();
            for (int i = 0; i < n_points; i++) {
                D0_full.block(2*i,3,1,3) = -x[i].homogeneous().transpose();
                D0_full.block(2*i,6,1,3) = y[i](1) * x[i].homogeneous().transpose();
                D0_full.block(2*i+1,0,1,3) = x[i].homogeneous().transpose();
                D0_full.block(2*i+1,6,1,3) = -y[i](0) * x[i].homogeneous().transpose();
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(D0_full, Eigen::ComputeFullV);
            Eigen::VectorXd h = svd.matrixV().col(8);
            Eigen::Matrix3d H = Eigen::Map<Eigen::Matrix3d>(h.data(), 3, 3).transpose();
            double res = 0.0;
            for (int j = 0; j < n_points; j++) {
                Eigen::Vector3d tmp = H * x[j].homogeneous();
                if (tmp(2) == 0.0) { res = std::numeric_limits<double>::max(); break; }
                res += (y[j] - tmp.hnormalized()).squaredNorm();
            }
            if (res < min_res) {
                min_res = res;
                posedata.homography = H;
                posedata.distortion_parameter2 = 0.0;
            }
        }

        if (min_res < std::numeric_limits<double>::max())
            output.push_back(posedata);

        return output;
    }
}  // namespace FitzgibbonCVPR2001
}  // namespace HomLib
