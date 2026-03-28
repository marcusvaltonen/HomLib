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
#include <float.h>  // For DBL_MAX
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include "posedata.hpp"
#include "radial.hpp"

namespace HomLib {
namespace FitzgibbonCVPR2001 {
    inline Eigen::Matrix3d vec2asym(const Eigen::Vector3d& t);

    std::vector<HomLib::PoseData> get(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y
    ) {
        // This is a five point method
        const int n_points = 5;

        // Make homogenous
        Eigen::Matrix<double, 3, n_points> x1, x2;
        for (int i = 0; i < n_points; i++) {
            x1.col(i) = x[i].homogeneous();
            x2.col(i) = y[i].homogeneous();
        }

        // Compute the distance to center point
        Eigen::Matrix<double, 3, n_points> z1;
        z1.setZero();
        Eigen::Matrix<double, 3, n_points> z2;
        z2.setZero();
        for (int i = 0; i < n_points; i++) {
            z1(2, i) = x[i].squaredNorm();
            z2(2, i) = y[i].squaredNorm();
        }

        // Initialize D0, D1 and D2
        Eigen::MatrixXd D0(9, 9);
        D0.setZero();
        Eigen::MatrixXd D1(9, 9);
        D1.setZero();
        Eigen::MatrixXd D2(9, 9);
        D2.setZero();

        Eigen::Matrix3d Bx2, Bz2;
        Eigen::Matrix3d e1, e2;
        Eigen::Vector3d tmp;

        for (int k = 0; k < n_points; k++) {
            tmp = x2.col(k);
            Bx2 = HomLib::FitzgibbonCVPR2001::vec2asym(tmp);
            tmp = z2.col(k);
            Bz2 = HomLib::FitzgibbonCVPR2001::vec2asym(tmp);

            // D0
            e1 = Bx2.row(0).transpose() * x1.col(k).transpose();
            e2 = Bx2.row(1).transpose() * x1.col(k).transpose();
            D0.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
            if (k < 4) {  // Assure it is 9x9
                D0.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
            }
            // D1
            e1 = Bx2.row(0).transpose() * z1.col(k).transpose() + Bz2.row(0).transpose() * x1.col(k).transpose();
            e2 = Bx2.row(1).transpose() * z1.col(k).transpose() + Bz2.row(1).transpose() * x1.col(k).transpose();
            D1.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
            if (k < 4) {
                D1.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
            }

            // D2
            e1 = Bz2.row(0).transpose() * z1.col(k).transpose();
            e2 = Bz2.row(1).transpose() * z1.col(k).transpose();
            D2.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
            if (k < 4) {
                D2.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
            }
        }

        // Create generalized eigenvalue problem
        Eigen::MatrixXd A(18, 18);
        Eigen::MatrixXd B(18, 18);
        A.setZero();
        B.setZero();

        A.topLeftCorner(9, 9) = -D0;
        A.bottomRightCorner(9, 9) = Eigen::MatrixXd::Identity(9, 9);
        B.topLeftCorner(9, 9) = D1;
        B.topRightCorner(9, 9) = D2;
        B.bottomLeftCorner(9, 9) = Eigen::MatrixXd::Identity(9, 9);

        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
        ges.compute(A, B, true);
        Eigen::VectorXcd l;
        Eigen::MatrixXd X;
        l = ges.eigenvalues();
        Eigen::MatrixXcd eigvecs;
        eigvecs = ges.eigenvectors();
        X = eigvecs.real().topRows(9);

        // Extract correct solution
        Eigen::Matrix3d Htmp;
        double minres = DBL_MAX;
        std::vector<HomLib::PoseData> output;
        HomLib::PoseData posedata;
        double ltmp;
        Eigen::Array<bool, 1, 18> is_ok;
        is_ok = l.array().isFinite() && l.array().imag() == 0;

        for (int k = 0; k < 18; k++) {
            if (is_ok(k)) {
                ltmp = l(k).real();
                Htmp = Eigen::Map<Eigen::Matrix3d>(X.col(k).data(), 3, 3);
                // Reproject points
                std::vector<Eigen::Vector2d> xu = radialundistort(x, ltmp);
                std::vector<Eigen::Vector2d> yu_est;
                for (int i = 0; i < n_points; i++) {
                    Eigen::Vector3d tmp = Htmp * xu[i].homogeneous();
                    yu_est.push_back(tmp.hnormalized());
                }
                std::vector<Eigen::Vector2d> y_est = radialdistort(yu_est, ltmp);

                // Compute residual in distorted space
                double res = 0.0;
                for (int i = 0; i < n_points; i++) {
                    res += (y[i] - y_est[i]).squaredNorm();
                }
                if (res < minres) {
                    minres = res;
                    posedata.homography = Htmp;
                    posedata.distortion_parameter = ltmp;
                }
            }
        }
        if (minres < DBL_MAX)
            output.push_back(posedata);

        return output;
    }

    // TODO(marcusvaltonen): Refactor -> helpers when necessary
    inline Eigen::Matrix3d vec2asym(const Eigen::Vector3d& t) {
        Eigen::Matrix3d t_hat;
        t_hat << 0, -t(2), t(1),
                 t(2), 0, -t(0),
                -t(1), t(0), 0;
        return t_hat;
    }
}  // namespace FitzgibbonCVPR2001
}  // namespace HomLib
