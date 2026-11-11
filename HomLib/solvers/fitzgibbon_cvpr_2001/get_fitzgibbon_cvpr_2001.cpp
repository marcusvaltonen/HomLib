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
        const std::vector<Eigen::Vector2d> &y,
        const HomLib::DistortionRange &range
    ) {
        const int n_points = static_cast<int>(x.size());
        if (n_points < 5)
            return {};

        // Build 9x9 QEP coefficient matrices C0, C1, C2.
        // Minimal (n_points == 5): fixed-size matrices, C = D directly.
        // Overdetermined (n_points > 5): QR-project the QEP onto D0 (avoids squaring its
        // condition number). Still biased with noise; refine the result.
        Eigen::Matrix<double, 9, 9> C0, C1, C2;

        if (n_points == 5) {
            Eigen::Matrix<double, 3, 5> x1, x2;
            for (int i = 0; i < 5; i++) {
                x1.col(i) = x[i].homogeneous();
                x2.col(i) = y[i].homogeneous();
            }
            Eigen::Matrix<double, 3, 5> z1, z2;
            z1.setZero();
            z2.setZero();
            for (int i = 0; i < 5; i++) {
                z1(2, i) = x[i].squaredNorm();
                z2(2, i) = y[i].squaredNorm();
            }

            Eigen::MatrixXd D0(9, 9), D1(9, 9), D2(9, 9);
            D0.setZero();
            D1.setZero();
            D2.setZero();
            Eigen::Matrix3d Bx2, Bz2, e1, e2;
            Eigen::Vector3d tmp;
            for (int k = 0; k < 5; k++) {
                tmp = x2.col(k);
                Bx2 = HomLib::FitzgibbonCVPR2001::vec2asym(tmp);
                tmp = z2.col(k);
                Bz2 = HomLib::FitzgibbonCVPR2001::vec2asym(tmp);
                e1 = Bx2.row(0).transpose() * x1.col(k).transpose();
                e2 = Bx2.row(1).transpose() * x1.col(k).transpose();
                D0.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
                if (k < 4) {
                    D0.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
                }
                e1 = Bx2.row(0).transpose() * z1.col(k).transpose() + Bz2.row(0).transpose() * x1.col(k).transpose();
                e2 = Bx2.row(1).transpose() * z1.col(k).transpose() + Bz2.row(1).transpose() * x1.col(k).transpose();
                D1.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
                if (k < 4) {
                    D1.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
                }
                e1 = Bz2.row(0).transpose() * z1.col(k).transpose();
                e2 = Bz2.row(1).transpose() * z1.col(k).transpose();
                D2.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
                if (k < 4) {
                    D2.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
                }
            }
            C0 = D0;
            C1 = D1;
            C2 = D2;
        } else {
            const int n_rows = 2 * n_points;
            Eigen::Matrix<double, 3, Eigen::Dynamic> x1(3, n_points), x2(3, n_points);
            for (int i = 0; i < n_points; i++) {
                x1.col(i) = x[i].homogeneous();
                x2.col(i) = y[i].homogeneous();
            }
            Eigen::Matrix<double, 3, Eigen::Dynamic> z1(3, n_points), z2(3, n_points);
            z1.setZero();
            z2.setZero();
            for (int i = 0; i < n_points; i++) {
                z1(2, i) = x[i].squaredNorm();
                z2(2, i) = y[i].squaredNorm();
            }

            Eigen::MatrixXd D0(n_rows, 9), D1(n_rows, 9), D2(n_rows, 9);
            D0.setZero();
            D1.setZero();
            D2.setZero();
            Eigen::Matrix3d Bx2, Bz2, e1, e2;
            Eigen::Vector3d tmp;
            for (int k = 0; k < n_points; k++) {
                tmp = x2.col(k);
                Bx2 = HomLib::FitzgibbonCVPR2001::vec2asym(tmp);
                tmp = z2.col(k);
                Bz2 = HomLib::FitzgibbonCVPR2001::vec2asym(tmp);
                e1 = Bx2.row(0).transpose() * x1.col(k).transpose();
                e2 = Bx2.row(1).transpose() * x1.col(k).transpose();
                D0.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
                D0.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
                e1 = Bx2.row(0).transpose() * z1.col(k).transpose() + Bz2.row(0).transpose() * x1.col(k).transpose();
                e2 = Bx2.row(1).transpose() * z1.col(k).transpose() + Bz2.row(1).transpose() * x1.col(k).transpose();
                D1.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
                D1.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
                e1 = Bz2.row(0).transpose() * z1.col(k).transpose();
                e2 = Bz2.row(1).transpose() * z1.col(k).transpose();
                D2.row(2*k) = Eigen::Map<Eigen::VectorXd>(e1.data(), 9);
                D2.row(2*k + 1) = Eigen::Map<Eigen::VectorXd>(e2.data(), 9);
            }
            // QR projection: avoids squaring D0's condition number (normal
            // equations fail when k ≈ 0 because D0 becomes rank-deficient).
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(D0);
            Eigen::MatrixXd thin_Q = qr.householderQ()
                * Eigen::MatrixXd::Identity(n_rows, 9);
            C0.setZero();
            C0.triangularView<Eigen::Upper>() = qr.matrixQR().topRows(9);
            C1 = thin_Q.transpose() * D1;
            C2 = thin_Q.transpose() * D2;
        }

        // Companion linearization of the QEP
        Eigen::MatrixXd A(18, 18);
        Eigen::MatrixXd B(18, 18);
        A.setZero();
        B.setZero();

        A.topLeftCorner(9, 9) = -C0;
        A.bottomRightCorner(9, 9) = Eigen::MatrixXd::Identity(9, 9);
        B.topLeftCorner(9, 9) = C1;
        B.topRightCorner(9, 9) = C2;
        B.bottomLeftCorner(9, 9) = Eigen::MatrixXd::Identity(9, 9);

        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
        ges.compute(A, B, true);
        if (ges.info() != Eigen::Success) {
            return {};
        }
        Eigen::VectorXcd l = ges.eigenvalues();
        Eigen::MatrixXd X = ges.eigenvectors().real().topRows(9);

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
                if (!range.contains(ltmp)) {
                    continue;
                }
                Htmp = Eigen::Map<Eigen::Matrix3d>(X.col(k).data(), 3, 3);
                // Reproject the points; a root for which a point has no reprojection (no
                // undistorted point, H x at infinity, no distorted point) is skipped
                double res = 0.0;
                bool ok = true;
                for (int i = 0; i < n_points && ok; i++) {
                    Eigen::Vector2d xu, y_est;
                    if (!HomLib::try_radialundistort(x[i], ltmp, xu)) {
                        ok = false;
                        break;
                    }
                    const Eigen::Vector3d tmp = Htmp * xu.homogeneous();
                    ok = tmp(2) != 0.0 && HomLib::try_radialdistort(tmp.hnormalized(), ltmp, y_est);
                    if (ok)
                        res += (y[i] - y_est).squaredNorm();   // residual in distorted space
                }
                if (!ok) {
                    continue;
                }
                if (res < minres) {
                    minres = res;
                    posedata.homography = Htmp;
                    posedata.distortion_parameter = ltmp;
                    posedata.distortion_parameter2 = ltmp;
                }
            }
        }
        // k = 0 fallback: the QEP is singular there, so recover the DLT
        // solution from the nullspace of D0.
        if (range.contains(0.0)) {
            Eigen::MatrixXd D0_full(2 * n_points, 9);
            D0_full.setZero();
            for (int i = 0; i < n_points; i++) {
                Eigen::RowVector3d xh = x[i].homogeneous().transpose();
                D0_full.block(2*i,   3, 1, 3) = -xh;
                D0_full.block(2*i,   6, 1, 3) = y[i](1) * xh;
                D0_full.block(2*i+1, 0, 1, 3) = xh;
                D0_full.block(2*i+1, 6, 1, 3) = -y[i](0) * xh;
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(D0_full, Eigen::ComputeFullV);
            Eigen::VectorXd h = svd.matrixV().col(8);
            Eigen::Matrix3d H_dlt = Eigen::Map<Eigen::Matrix3d>(h.data(), 3, 3).transpose();
            double res = 0.0;
            bool ok = true;
            for (int i = 0; i < n_points && ok; i++) {
                const Eigen::Vector3d tmp = H_dlt * x[i].homogeneous();
                ok = tmp(2) != 0.0;
                if (ok)
                    res += (y[i] - tmp.hnormalized()).squaredNorm();
            }
            if (ok && res < minres) {
                minres = res;
                posedata.homography = H_dlt;
                posedata.distortion_parameter = 0.0;
                posedata.distortion_parameter2 = 0.0;
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
