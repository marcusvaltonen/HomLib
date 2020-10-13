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

#include "get_valtonenornhag_arxiv_2020b.hpp"
#include <Eigen/Geometry>
#include <cmath>  // max
#include <vector>
#include "solver_valtonenornhag_arxiv_2020b_fHf.hpp"
#include "normalize2dpts.hpp"
#include "posedata.hpp"

namespace HomLib::ValtonenOrnhagArxiv2020B {
    inline Eigen::Vector4d construct_hvector(double w, const Eigen::VectorXd input);

    std::vector<HomLib::PoseData> get_fHf(
        const Eigen::MatrixXd &p1,
        const Eigen::MatrixXd &p2,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &R2
    ) {
        // This is a 2-point method
        int nbr_pts = 2;

        // We expect inhomogenous input data, i.e. p1 and p2 are 2x3 matrices
        assert(p1.rows() == 2);
        assert(p2.rows() == 2);
        assert(p1.cols() == nbr_pts);
        assert(p2.cols() == nbr_pts);
        int nbr_coeffs = 26;

        // Save copies of the inverse rotation
        Eigen::Matrix3d R1T = R1.transpose();
        Eigen::Matrix3d R2T = R2.transpose();

        // Compute normalization matrix
        double scale1 = normalize2dpts(p1);
        double scale2 = normalize2dpts(p2);
        double scale = std::max(scale1, scale2);
        Eigen::Vector3d s;
        s << scale, scale, 1.0;
        Eigen::DiagonalMatrix<double, 3> S = s.asDiagonal();

        // Normalize data
        Eigen::MatrixXd x1(3, 2);
        Eigen::MatrixXd x2(3, 2);
        x1 = p1.colwise().homogeneous();
        x2 = p2.colwise().homogeneous();

        x1 = S * x1;
        x2 = S * x2;

        Eigen::Matrix2d x1t;
        x1t << x1.colwise().hnormalized();
        Eigen::Matrix2d x2t;
        x2t << x2.colwise().hnormalized();

        // Wrap input data to expected format
        Eigen::VectorXd input(nbr_coeffs);
        input << x1t.col(0),
                 x2t.col(0),
                 x1t.col(1),
                 x2t.col(1),
                 Eigen::Map<Eigen::VectorXd>(R1T.data(), 9),
                 Eigen::Map<Eigen::VectorXd>(R2T.data(), 9);

        // Extract w
        Eigen::VectorXcd w = HomLib::ValtonenOrnhagArxiv2020B::solver_fHf(input);

        // Pre-processing: Remove complex-valued solutions
        double thresh = 1e-5;
        Eigen::ArrayXd real_w = w.imag().array().abs();

        // This is a 2 pt solver
        std::vector<HomLib::PoseData> posedata;
        HomLib::PoseData tmp_pose;
        double w_tmp;
        Eigen::Vector4d hvec;
        Eigen::Matrix3d K, Ki, Htmp;

        double tol = 1e-13;

        for (int i = 0; i < real_w.size(); i++) {
            if (real_w(i) <= thresh) {
                // Compute algebraic error, and compare to other solutions.
                w_tmp = w(i).real();

                // Compute h vector
                hvec = construct_hvector(w_tmp, input);

                // Two spurious solutions were added, corresponding to last element
                // equal to zero.
                if (std::abs(hvec(3)) > tol) {
                    Htmp = Eigen::Matrix3d::Identity(3, 3);
                    Htmp.col(1) = hvec.hnormalized();

                    K = Eigen::Vector3d(w_tmp, w_tmp, 1).asDiagonal();
                    Ki = Eigen::Vector3d(1, 1, w_tmp).asDiagonal();
                    Htmp = S.inverse() * K * R2 * Htmp * R1T * Ki * S;

                    // Package output
                    tmp_pose.homography = Htmp;
                    tmp_pose.focal_length = w_tmp / scale;
                    posedata.push_back(tmp_pose);
                }
            }
        }

        return posedata;
    }

    inline Eigen::Vector4d construct_hvector(double w, const Eigen::VectorXd input) {
        Eigen::VectorXd d(27);
        d << w, input;

        double s1 = d[0]*d[25];
        double s2 = d[0]*d[16];
        double t1 = s1 + d[7]*d[19] + d[8]*d[22];
        double t2 = d[3]*d[19] + s1 + d[4]*d[22];
        double t3 = d[5]*d[10] + s2 + d[6]*d[13];
        double t4 = d[1]*d[10] + d[2]*d[13] + s2;
        double u1 = d[0]*d[24];
        double u2 = d[0]*d[15];
        double u3 = d[0]*d[26];

        Eigen::Matrix<double, 4, 4> M;
        M << 0, -t4*(d[3]*d[20] + u3 + d[4]*d[23]), t4*t2, (d[1]*d[11] + d[2]*d[14] + d[0]*d[17])*t2,  // NOLINT
             -t4*t2,  t4*(d[3]*d[18] + u1 + d[4]*d[21]), 0, -(d[1]*d[9] + d[2]*d[12] + u2)*t2,  // NOLINT
             0, -t3*(u3 + d[7]*d[20] + d[8]*d[23]), t3*t1,  (d[5]*d[11] + d[0]*d[17] + d[6]*d[14])*t1,  // NOLINT
             -t3*t1,  t3*(u1 + d[7]*d[18] + d[8]*d[21]), 0, -(d[5]*d[9] + u2 + d[6]*d[12])*t1;  // NOLINT

        // Perform SVD
        Eigen::JacobiSVD<Eigen::Matrix<double, 4, 4>> svd(M, Eigen::ComputeFullV);

        // Extract hvector
        Eigen::Vector4d h = svd.matrixV().col(3);

        return h;
    }
}  // namespace HomLib::ValtonenOrnhagArxiv2020B
