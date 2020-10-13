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

#include "get_valtonenornhag_arxiv_2020a.hpp"
#include <float.h>  // DBL_MAX
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>  // abs
#include "solver_valtonenornhag_arxiv_2020a_fHf.hpp"
#include "normalize2dpts.hpp"

namespace HomLib::ValtonenOrnhagArxiv2020A {

    inline double get_algebraic_error_floor_fHf(const Eigen::VectorXd &data);

    HomLib::PoseData get_fHf(
        const Eigen::MatrixXd &p1,
        const Eigen::MatrixXd &p2,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &R2
    ) {
        int nbr_coeffs = 30;
        int nbr_unknowns = 6;

        // Save copies of the inverse rotation
        Eigen::Matrix3d R1T = R1.transpose();
        Eigen::Matrix3d R2T = R2.transpose();

        // Compute normalization matrix
        double scale = normalize2dpts(p1);
        Eigen::Vector3d s;
        s << scale, scale, 1.0;
        Eigen::DiagonalMatrix<double, 3> S = s.asDiagonal();

        // Normalize data
        Eigen::Matrix3d x1;
        Eigen::Matrix3d x2;
        x1 = p1.colwise().homogeneous();
        x2 = p2.colwise().homogeneous();

        x1 = S * x1;
        x2 = S * x2;

        Eigen::MatrixXd x1t(2, 3);
        x1t << x1.colwise().hnormalized();
        Eigen::MatrixXd x2t(2, 3);
        x2t << x2.colwise().hnormalized();

        // Wrap input data to expected format
        Eigen::VectorXd input(nbr_coeffs);
        input << x1t.col(0),
                 x2t.col(0),
                 x1t.col(1),
                 x2t.col(1),
                 x1t.col(2),
                 x2t.col(2),
                 Eigen::Map<Eigen::VectorXd>(R1T.data(), 9),
                 Eigen::Map<Eigen::VectorXd>(R2T.data(), 9);

        // Extract solution
        Eigen::MatrixXcd sols = HomLib::ValtonenOrnhagArxiv2020A::solver_fHf(input);

        // Pre-processing: Remove complex-valued solutions
        double thresh = 1e-5;
        Eigen::ArrayXd real_sols(7);
        real_sols = sols.imag().cwiseAbs().colwise().sum();
        int nbr_real_sols = (real_sols <= thresh).count();

        // Allocate space for putative (real) homographies
        Eigen::MatrixXd best_homography(3, 3);
        double best_focal_length;
        double best_algebraic_error = DBL_MAX;
        double algebraic_error;

        // Since this is a 2.5 pt solver, use the last
        // (previously unused) constraint, to discard
        // false solutions.
        Eigen::ArrayXd xx(6);
        Eigen::VectorXd input_algebraic(nbr_coeffs + nbr_unknowns);

        for (int i = 0; i < real_sols.size(); i++) {
            if (real_sols(i) <= thresh) {
                // Compute algebraic error, and compare to other solutions.
                xx = sols.col(i).real();
                input_algebraic << xx, input;
                algebraic_error = HomLib::ValtonenOrnhagArxiv2020A::get_algebraic_error_floor_fHf(input_algebraic);

                if (algebraic_error < best_algebraic_error) {
                    best_algebraic_error = algebraic_error;
                    best_homography << xx[0], xx[2], xx[1],
                                           0, xx[3],     0,
                                      -xx[1], xx[4], xx[0];
                    best_focal_length = xx[5];
                }
            }
        }
        // Construct homography
        Eigen::Matrix3d K, Ki, H;
        K = Eigen::Vector3d(best_focal_length, best_focal_length, 1).asDiagonal();
        Ki = Eigen::Vector3d(1, 1, best_focal_length).asDiagonal();
        H = S.inverse() * K * R2 * best_homography * R1.transpose() * Ki * S;

        // Package output
        HomLib::PoseData posedata;
        posedata.homography = H;
        posedata.focal_length = best_focal_length / scale;

        return posedata;
    }

    // Function that utilizes the last equation of the DLT system to discard false solutions
    inline double get_algebraic_error_floor_fHf(const Eigen::VectorXd &data) {
        const double* d = data.data();

        // Compute algebraic error
        double error;
        error = -d[0]*std::pow(d[5], 2)*d[24]*d[34] - d[0]*d[5]*d[14]*d[18]*d[34] - d[0]*d[5]*d[15]*d[21]*d[34]
            - d[0]*d[5]*d[16]*d[24]*d[28] - d[0]*d[5]*d[17]*d[24]*d[31] - d[0]*d[14]*d[16]*d[18]*d[28]
            - d[0]*d[14]*d[17]*d[18]*d[31] - d[0]*d[15]*d[16]*d[21]*d[28] - d[0]*d[15]*d[17]*d[21]*d[31]
            - d[1]*std::pow(d[5], 2)*d[26]*d[34] - d[1]*d[5]*d[14]*d[20]*d[34] - d[1]*d[5]*d[15]*d[23]*d[34]
            - d[1]*d[5]*d[16]*d[26]*d[28] - d[1]*d[5]*d[17]*d[26]*d[31] - d[1]*d[14]*d[16]*d[20]*d[28]
            - d[1]*d[14]*d[17]*d[20]*d[31] - d[1]*d[15]*d[16]*d[23]*d[28] - d[1]*d[15]*d[17]*d[23]*d[31]
            - d[2]*std::pow(d[5], 2)*d[25]*d[34] - d[2]*d[5]*d[14]*d[19]*d[34] - d[2]*d[5]*d[15]*d[22]*d[34]
            - d[2]*d[5]*d[16]*d[25]*d[28] - d[2]*d[5]*d[17]*d[25]*d[31] - d[2]*d[14]*d[16]*d[19]*d[28]
            - d[2]*d[14]*d[17]*d[19]*d[31] - d[2]*d[15]*d[16]*d[22]*d[28] - d[2]*d[15]*d[17]*d[22]*d[31]
            + d[3]*std::pow(d[5], 2)*d[25]*d[33] + d[3]*d[5]*d[14]*d[19]*d[33] + d[3]*d[5]*d[15]*d[22]*d[33]
            + d[3]*d[5]*d[16]*d[25]*d[27] + d[3]*d[5]*d[17]*d[25]*d[30] + d[3]*d[14]*d[16]*d[19]*d[27]
            + d[3]*d[14]*d[17]*d[19]*d[30] + d[3]*d[15]*d[16]*d[22]*d[27] + d[3]*d[15]*d[17]*d[22]*d[30];
        return abs(error);
    }
}  // namespace HomLib::ValtonenOrnhagArxiv2020A
