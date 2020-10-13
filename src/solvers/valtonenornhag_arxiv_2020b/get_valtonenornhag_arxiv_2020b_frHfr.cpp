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
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <float.h>  // DBL_MAX
#include <vector>
#include <cmath>  // abs
#include "solver_valtonenornhag_arxiv_2020b_frHfr.hpp"
#include "normalize2dpts.hpp"
#include "posedata.hpp"
#include "gj.hpp"


namespace HomLib::ValtonenOrnhagArxiv2020B {
    inline Eigen::Matrix<double, 5, 1> construct_sols(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& input,
        const Eigen::MatrixXd& M);
    inline Eigen::Vector4d rot2quat(const Eigen::Matrix3d& R);
    inline double get_algebraic_error_norot_frHfr(const Eigen::VectorXd &data);

    HomLib::PoseData get_frHfr(
        const Eigen::MatrixXd &p1,
        const Eigen::MatrixXd &p2,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d & R2
    ) {
        // This is a 2.5 point method
        const int nbr_pts = 3;

        // We expect inhomogenous input data, i.e. p1 and p2 are 2x5 matrices
        assert(p1.rows() == 2);
        assert(p2.rows() == 2);
        assert(p1.cols() == nbr_pts);
        assert(p2.cols() == nbr_pts);

        // Compute normalization matrix
        double scale1 = HomLib::normalize2dpts(p1);
        double scale2 = HomLib::normalize2dpts(p2);
        double scale = std::max(scale1, scale2);
        Eigen::Vector3d s;
        s << scale, scale, 1.0;
        Eigen::DiagonalMatrix<double, 3> S = s.asDiagonal();

        // Normalize data
        Eigen::MatrixXd x1(3, nbr_pts);
        Eigen::MatrixXd x2(3, nbr_pts);
        x1 = p1.colwise().homogeneous();
        x2 = p2.colwise().homogeneous();

        x1 = S * x1;
        x2 = S * x2;

        Eigen::MatrixXd u1(2, nbr_pts);
        u1 << x1.colwise().hnormalized();
        Eigen::MatrixXd u2(2, nbr_pts);
        u2 << x2.colwise().hnormalized();

        // Save copies of the modified rotation matrix
        Eigen::Matrix3d R = R2 * R1.transpose();
        Eigen::Vector4d q = rot2quat(R);

        // Wrap input data to expected format
        Eigen::VectorXd d(19);
        d << u1.col(0),
             u2.col(0),
             u1.col(1),
             u2.col(1),
             u1.col(2),
             u2.col(2),
             q,
             R1.col(1);

        // Create M matrix
        // TODO(marcusvaltonen): Optimize this
        Eigen::MatrixXd M(3, 9);
        M << std::pow(d[0],2)*d[2]*d[18] + std::pow(d[1],2)*d[2]*d[18], d[2]*d[18], d[0]*d[2]*d[16] + d[1]*d[2]*d[17], -std::pow(d[0],2)*d[3]*d[18] - std::pow(d[1],2)*d[3]*d[18], -d[3]*d[18], -d[0]*d[3]*d[16] - d[1]*d[3]*d[17], -2*std::pow(d[0],2)*d[2]*d[12]*d[13] + 2*std::pow(d[0],2)*d[2]*d[14]*d[15] - 2*std::pow(d[0],2)*d[3]*d[12]*d[14] - 2*std::pow(d[0],2)*d[3]*d[13]*d[15] - 2*std::pow(d[1],2)*d[2]*d[12]*d[13] + 2*std::pow(d[1],2)*d[2]*d[14]*d[15] - 2*std::pow(d[1],2)*d[3]*d[12]*d[14] - 2*std::pow(d[1],2)*d[3]*d[13]*d[15], -2*d[2]*d[12]*d[13] + 2*d[2]*d[14]*d[15] - 2*d[3]*d[12]*d[14] - 2*d[3]*d[13]*d[15], 2*d[0]*d[2]*d[12]*d[15] + 2*d[0]*d[2]*d[13]*d[14] + 2*d[0]*d[3]*std::pow(d[14],2) + 2*d[0]*d[3]*std::pow(d[15],2) - d[0]*d[3] - 2*d[1]*d[2]*std::pow(d[13],2) - 2*d[1]*d[2]*std::pow(d[15],2) + d[1]*d[2] + 2*d[1]*d[3]*d[12]*d[15] - 2*d[1]*d[3]*d[13]*d[14],  // NOLINT
        std::pow(d[4],2)*d[6]*d[18] + std::pow(d[5],2)*d[6]*d[18], d[6]*d[18], d[4]*d[6]*d[16] + d[5]*d[6]*d[17], -std::pow(d[4],2)*d[7]*d[18] - std::pow(d[5],2)*d[7]*d[18], -d[7]*d[18], -d[4]*d[7]*d[16] - d[5]*d[7]*d[17], -2*std::pow(d[4],2)*d[6]*d[12]*d[13] + 2*std::pow(d[4],2)*d[6]*d[14]*d[15] - 2*std::pow(d[4],2)*d[7]*d[12]*d[14] - 2*std::pow(d[4],2)*d[7]*d[13]*d[15] - 2*std::pow(d[5],2)*d[6]*d[12]*d[13] + 2*std::pow(d[5],2)*d[6]*d[14]*d[15] - 2*std::pow(d[5],2)*d[7]*d[12]*d[14] - 2*std::pow(d[5],2)*d[7]*d[13]*d[15], -2*d[6]*d[12]*d[13] + 2*d[6]*d[14]*d[15] - 2*d[7]*d[12]*d[14] - 2*d[7]*d[13]*d[15], 2*d[4]*d[6]*d[12]*d[15] + 2*d[4]*d[6]*d[13]*d[14] + 2*d[4]*d[7]*std::pow(d[14],2) + 2*d[4]*d[7]*std::pow(d[15],2) - d[4]*d[7] - 2*d[5]*d[6]*std::pow(d[13],2) - 2*d[5]*d[6]*std::pow(d[15],2) + d[5]*d[6] + 2*d[5]*d[7]*d[12]*d[15] - 2*d[5]*d[7]*d[13]*d[14],  // NOLINT
        std::pow(d[8],2)*d[10]*d[18] + std::pow(d[9],2)*d[10]*d[18], d[10]*d[18], d[8]*d[10]*d[16] + d[9]*d[10]*d[17], -std::pow(d[8],2)*d[11]*d[18] - std::pow(d[9],2)*d[11]*d[18], -d[11]*d[18], -d[8]*d[11]*d[16] - d[9]*d[11]*d[17], -2*std::pow(d[8],2)*d[10]*d[12]*d[13] + 2*std::pow(d[8],2)*d[10]*d[14]*d[15] - 2*std::pow(d[8],2)*d[11]*d[12]*d[14] - 2*std::pow(d[8],2)*d[11]*d[13]*d[15] - 2*std::pow(d[9],2)*d[10]*d[12]*d[13] + 2*std::pow(d[9],2)*d[10]*d[14]*d[15] - 2*std::pow(d[9],2)*d[11]*d[12]*d[14] - 2*std::pow(d[9],2)*d[11]*d[13]*d[15], -2*d[10]*d[12]*d[13] + 2*d[10]*d[14]*d[15] - 2*d[11]*d[12]*d[14] - 2*d[11]*d[13]*d[15], 2*d[8]*d[10]*d[12]*d[15] + 2*d[8]*d[10]*d[13]*d[14] + 2*d[8]*d[11]*std::pow(d[14],2) + 2*d[8]*d[11]*std::pow(d[15],2) - d[8]*d[11] - 2*d[9]*d[10]*std::pow(d[13],2) - 2*d[9]*d[10]*std::pow(d[15],2) + d[9]*d[10] + 2*d[9]*d[11]*d[12]*d[15] - 2*d[9]*d[11]*d[13]*d[14];  // NOLINT

        // Find row-reduced echelon form
        HomLib::gj(&M);

        // Wrap input data to expected format
        Eigen::VectorXd input(37);
        input << Eigen::Map<Eigen::VectorXd>(M.rightCols(6).data(), 6*3), d;

        // Extract solution
        Eigen::MatrixXcd sols = HomLib::ValtonenOrnhagArxiv2020B::solver_frHfr(input);

        // Pre-processing: Remove complex-valued solutions
        double thresh = 1e-5;
        Eigen::Array3d real_sols;
        real_sols = sols.imag().cwiseAbs().colwise().sum();

        // Allocate space for putative (real) homographies
        Eigen::Matrix<double, 5, 1> this_sols;
        Eigen::Matrix<double, 5, 1> best_sols;
        double best_algebraic_error = DBL_MAX;
        double algebraic_error;

        // Since this is a 2.5 pt solver, use the last
        // (previously unused) constraint, to discard
        // false solutions.
        Eigen::ArrayXd xx(2);
        Eigen::VectorXd input_algebraic(29);
        input_algebraic << Eigen::VectorXd::Zero(5),
                           u1.col(0),
                           u2.col(0),
                           u1.col(1),
                           u2.col(1),
                           u1.col(2),
                           u2.col(2),
                           Eigen::Map<Eigen::VectorXd>(R.data(), 9),
                           R1.col(1);

        for (int i = 0; i < real_sols.size(); i++) {
            if (real_sols(i) <= thresh) {
                // Get parameters.
                xx = sols.col(i).real();

                // Construct complete set of motion parameters and intrinsic
                this_sols = construct_sols(xx, input, M.rightCols(6));

                // Test algebraic error
                input_algebraic.head(5) = this_sols;
                algebraic_error = get_algebraic_error_norot_frHfr(input_algebraic);

                if (algebraic_error < best_algebraic_error) {
                    best_algebraic_error = algebraic_error;
                    best_sols = this_sols;
                }
            }
        }

        // Construct output
        Eigen::Matrix3d K, Ki, H, Hy;
        Eigen::RowVector3d n(0, 1, 0);
        Hy = Eigen::Matrix3d::Identity(3, 3) + (R2.transpose() * best_sols.head(3)) * n;
        double f = best_sols(3);
        double r = best_sols(4);

        K = Eigen::Vector3d(f, f, 1).asDiagonal();
        Ki = Eigen::Vector3d(1, 1, f).asDiagonal();
        H = S.inverse() * K * R2 * Hy * R1.transpose() * Ki * S;

        // Account for scale
        HomLib::PoseData posedata;
        posedata.homography = H;
        posedata.focal_length = f / scale;
        posedata.distortion_parameter = r * std::pow(scale, 2);

        return posedata;
    }

    // Function that utilizes the last equation of the DLT system to discard false solutions
    inline Eigen::Matrix<double, 5, 1> construct_sols(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& input,
        const Eigen::MatrixXd& M
    ) {
        Eigen::VectorXd d(41);
        d << 0, 0, xx(2), xx(3), input;

        Eigen::Matrix<double, 5, 3> M2;
        // TODO(marcusvaltonen): Optimize this
        M2 << std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[22],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[22],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[23],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[23],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[22],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[23],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[22],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[22],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[23],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[23],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*d[3]*d[6]*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[22],2)*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[23],2)*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*d[9]*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[22]*std::pow(d[24],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[22]*std::pow(d[25],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[23]*std::pow(d[24],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[23]*std::pow(d[25],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[22],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[22],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[23],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[23],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],2)*d[3]*d[6]*d[22]*d[38] + std::pow(d[2],2)*d[3]*d[6]*d[23]*d[39] + std::pow(d[2],2)*d[3]*d[9]*d[22]*std::pow(d[24],2)*d[38] + std::pow(d[2],2)*d[3]*d[9]*d[22]*std::pow(d[25],2)*d[38] + std::pow(d[2],2)*d[3]*d[9]*d[23]*std::pow(d[24],2)*d[39] + std::pow(d[2],2)*d[3]*d[9]*d[23]*std::pow(d[25],2)*d[39] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[22],2)*d[40] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[23],2)*d[40] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[24],2)*d[40] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[25],2)*d[40] + std::pow(d[2],2)*d[9]*d[22]*d[38] + std::pow(d[2],2)*d[9]*d[23]*d[39] + std::pow(d[2],2)*d[12]*d[40] + d[2]*d[3]*d[12]*d[22]*std::pow(d[24],2)*d[38] + d[2]*d[3]*d[12]*d[22]*std::pow(d[25],2)*d[38] + d[2]*d[3]*d[12]*d[23]*std::pow(d[24],2)*d[39] + d[2]*d[3]*d[12]*d[23]*std::pow(d[25],2)*d[39] + d[2]*d[12]*d[22]*d[38] + d[2]*d[12]*d[23]*d[39], d[2]*d[3]*std::pow(d[22],2)*d[25]*d[40] + d[2]*d[3]*std::pow(d[23],2)*d[25]*d[40] + d[2]*d[25]*d[40] + d[22]*d[25]*d[38] + d[23]*d[25]*d[39], std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[22],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[22],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[23],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[23],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[22],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[23],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[22],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[22],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[23],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[23],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*d[3]*d[15]*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[22],2)*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[23],2)*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[24],2)*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[25],2)*d[40] + std::pow(d[2],3)*d[18]*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[22]*std::pow(d[24],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[22]*std::pow(d[25],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[23]*std::pow(d[24],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[23]*std::pow(d[25],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[22],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[22],2)*std::pow(d[25],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[23],2)*std::pow(d[24],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[23],2)*std::pow(d[25],2)*d[40] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[22],2)*std::pow(d[24],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[22],2)*std::pow(d[24],2)*d[36]*d[37] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[22],2)*std::pow(d[25],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[22],2)*std::pow(d[25],2)*d[36]*d[37] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[23],2)*std::pow(d[24],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[23],2)*std::pow(d[24],2)*d[36]*d[37] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[23],2)*std::pow(d[25],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[23],2)*std::pow(d[25],2)*d[36]*d[37] + std::pow(d[2],2)*d[3]*d[15]*d[22]*d[38] + std::pow(d[2],2)*d[3]*d[15]*d[23]*d[39] + std::pow(d[2],2)*d[3]*d[18]*d[22]*std::pow(d[24],2)*d[38] + std::pow(d[2],2)*d[3]*d[18]*d[22]*std::pow(d[25],2)*d[38] + std::pow(d[2],2)*d[3]*d[18]*d[23]*std::pow(d[24],2)*d[39] + std::pow(d[2],2)*d[3]*d[18]*d[23]*std::pow(d[25],2)*d[39] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[22],2)*d[40] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[23],2)*d[40] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[24],2)*d[40] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[25],2)*d[40] + 2*std::pow(d[2],2)*d[3]*std::pow(d[22],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[22],2)*d[36]*d[37] + 2*std::pow(d[2],2)*d[3]*std::pow(d[23],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[23],2)*d[36]*d[37] + 2*std::pow(d[2],2)*d[3]*std::pow(d[24],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[24],2)*d[36]*d[37] + 2*std::pow(d[2],2)*d[3]*std::pow(d[25],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[25],2)*d[36]*d[37] + std::pow(d[2],2)*d[18]*d[22]*d[38] + std::pow(d[2],2)*d[18]*d[23]*d[39] + std::pow(d[2],2)*d[21]*d[40] + 2*std::pow(d[2],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[36]*d[37] + d[2]*d[3]*d[21]*d[22]*std::pow(d[24],2)*d[38] + d[2]*d[3]*d[21]*d[22]*std::pow(d[25],2)*d[38] + d[2]*d[3]*d[21]*d[23]*std::pow(d[24],2)*d[39] + d[2]*d[3]*d[21]*d[23]*std::pow(d[25],2)*d[39] + d[2]*d[3]*std::pow(d[22],2)*d[25]*std::pow(d[34],2) - d[2]*d[3]*std::pow(d[22],2)*d[25]*std::pow(d[35],2) - d[2]*d[3]*std::pow(d[22],2)*d[25]*std::pow(d[36],2) + d[2]*d[3]*std::pow(d[22],2)*d[25]*std::pow(d[37],2) - 2*d[2]*d[3]*d[22]*std::pow(d[24],2)*d[34]*d[37] - 2*d[2]*d[3]*d[22]*std::pow(d[24],2)*d[35]*d[36] - 2*d[2]*d[3]*d[22]*std::pow(d[25],2)*d[34]*d[37] - 2*d[2]*d[3]*d[22]*std::pow(d[25],2)*d[35]*d[36] + d[2]*d[3]*std::pow(d[23],2)*d[25]*std::pow(d[34],2) - d[2]*d[3]*std::pow(d[23],2)*d[25]*std::pow(d[35],2) - d[2]*d[3]*std::pow(d[23],2)*d[25]*std::pow(d[36],2) + d[2]*d[3]*std::pow(d[23],2)*d[25]*std::pow(d[37],2) - d[2]*d[3]*d[23]*std::pow(d[24],2)*std::pow(d[34],2) + d[2]*d[3]*d[23]*std::pow(d[24],2)*std::pow(d[35],2) - d[2]*d[3]*d[23]*std::pow(d[24],2)*std::pow(d[36],2) + d[2]*d[3]*d[23]*std::pow(d[24],2)*std::pow(d[37],2) - d[2]*d[3]*d[23]*std::pow(d[25],2)*std::pow(d[34],2) + d[2]*d[3]*d[23]*std::pow(d[25],2)*std::pow(d[35],2) - d[2]*d[3]*d[23]*std::pow(d[25],2)*std::pow(d[36],2) + d[2]*d[3]*d[23]*std::pow(d[25],2)*std::pow(d[37],2) + d[2]*d[21]*d[22]*d[38] + d[2]*d[21]*d[23]*d[39] - 2*d[2]*d[22]*d[34]*d[37] - 2*d[2]*d[22]*d[35]*d[36] - d[2]*d[23]*std::pow(d[34],2) + d[2]*d[23]*std::pow(d[35],2) - d[2]*d[23]*std::pow(d[36],2) + d[2]*d[23]*std::pow(d[37],2) + d[2]*d[25]*std::pow(d[34],2) - d[2]*d[25]*std::pow(d[35],2) - d[2]*d[25]*std::pow(d[36],2) + d[2]*d[25]*std::pow(d[37],2) - 2*d[22]*d[25]*d[34]*d[36] + 2*d[22]*d[25]*d[35]*d[37] + 2*d[23]*d[25]*d[34]*d[35] + 2*d[23]*d[25]*d[36]*d[37],  // NOLINT
        std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[26],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[26],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[27],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[6]*std::pow(d[27],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[26],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[27],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[6]*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[26],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[26],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[27],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[9]*std::pow(d[27],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*d[3]*d[6]*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[26],2)*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[27],2)*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*d[3]*d[9]*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*d[9]*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[26]*std::pow(d[28],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[26]*std::pow(d[29],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[27]*std::pow(d[28],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[6]*d[27]*std::pow(d[29],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[26],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[26],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[27],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[12]*std::pow(d[27],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],2)*d[3]*d[6]*d[26]*d[38] + std::pow(d[2],2)*d[3]*d[6]*d[27]*d[39] + std::pow(d[2],2)*d[3]*d[9]*d[26]*std::pow(d[28],2)*d[38] + std::pow(d[2],2)*d[3]*d[9]*d[26]*std::pow(d[29],2)*d[38] + std::pow(d[2],2)*d[3]*d[9]*d[27]*std::pow(d[28],2)*d[39] + std::pow(d[2],2)*d[3]*d[9]*d[27]*std::pow(d[29],2)*d[39] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[26],2)*d[40] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[27],2)*d[40] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[28],2)*d[40] + std::pow(d[2],2)*d[3]*d[12]*std::pow(d[29],2)*d[40] + std::pow(d[2],2)*d[9]*d[26]*d[38] + std::pow(d[2],2)*d[9]*d[27]*d[39] + std::pow(d[2],2)*d[12]*d[40] + d[2]*d[3]*d[12]*d[26]*std::pow(d[28],2)*d[38] + d[2]*d[3]*d[12]*d[26]*std::pow(d[29],2)*d[38] + d[2]*d[3]*d[12]*d[27]*std::pow(d[28],2)*d[39] + d[2]*d[3]*d[12]*d[27]*std::pow(d[29],2)*d[39] + d[2]*d[12]*d[26]*d[38] + d[2]*d[12]*d[27]*d[39], d[2]*d[3]*std::pow(d[26],2)*d[29]*d[40] + d[2]*d[3]*std::pow(d[27],2)*d[29]*d[40] + d[2]*d[29]*d[40] + d[26]*d[29]*d[38] + d[27]*d[29]*d[39], std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[26],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[26],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[27],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],3)*d[15]*std::pow(d[27],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[26],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[27],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[15]*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[26],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[26],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[27],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*std::pow(d[3],2)*d[18]*std::pow(d[27],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*d[3]*d[15]*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[26],2)*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[27],2)*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[28],2)*d[40] + std::pow(d[2],3)*d[3]*d[18]*std::pow(d[29],2)*d[40] + std::pow(d[2],3)*d[18]*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[26]*std::pow(d[28],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[26]*std::pow(d[29],2)*d[38] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[27]*std::pow(d[28],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[15]*d[27]*std::pow(d[29],2)*d[39] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[26],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[26],2)*std::pow(d[29],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[27],2)*std::pow(d[28],2)*d[40] + std::pow(d[2],2)*std::pow(d[3],2)*d[21]*std::pow(d[27],2)*std::pow(d[29],2)*d[40] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[26],2)*std::pow(d[28],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[26],2)*std::pow(d[28],2)*d[36]*d[37] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[26],2)*std::pow(d[29],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[26],2)*std::pow(d[29],2)*d[36]*d[37] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[27],2)*std::pow(d[28],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[27],2)*std::pow(d[28],2)*d[36]*d[37] + 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[27],2)*std::pow(d[29],2)*d[34]*d[35] - 2*std::pow(d[2],2)*std::pow(d[3],2)*std::pow(d[27],2)*std::pow(d[29],2)*d[36]*d[37] + std::pow(d[2],2)*d[3]*d[15]*d[26]*d[38] + std::pow(d[2],2)*d[3]*d[15]*d[27]*d[39] + std::pow(d[2],2)*d[3]*d[18]*d[26]*std::pow(d[28],2)*d[38] + std::pow(d[2],2)*d[3]*d[18]*d[26]*std::pow(d[29],2)*d[38] + std::pow(d[2],2)*d[3]*d[18]*d[27]*std::pow(d[28],2)*d[39] + std::pow(d[2],2)*d[3]*d[18]*d[27]*std::pow(d[29],2)*d[39] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[26],2)*d[40] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[27],2)*d[40] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[28],2)*d[40] + std::pow(d[2],2)*d[3]*d[21]*std::pow(d[29],2)*d[40] + 2*std::pow(d[2],2)*d[3]*std::pow(d[26],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[26],2)*d[36]*d[37] + 2*std::pow(d[2],2)*d[3]*std::pow(d[27],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[27],2)*d[36]*d[37] + 2*std::pow(d[2],2)*d[3]*std::pow(d[28],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[28],2)*d[36]*d[37] + 2*std::pow(d[2],2)*d[3]*std::pow(d[29],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[3]*std::pow(d[29],2)*d[36]*d[37] + std::pow(d[2],2)*d[18]*d[26]*d[38] + std::pow(d[2],2)*d[18]*d[27]*d[39] + std::pow(d[2],2)*d[21]*d[40] + 2*std::pow(d[2],2)*d[34]*d[35] - 2*std::pow(d[2],2)*d[36]*d[37] + d[2]*d[3]*d[21]*d[26]*std::pow(d[28],2)*d[38] + d[2]*d[3]*d[21]*d[26]*std::pow(d[29],2)*d[38] + d[2]*d[3]*d[21]*d[27]*std::pow(d[28],2)*d[39] + d[2]*d[3]*d[21]*d[27]*std::pow(d[29],2)*d[39] + d[2]*d[3]*std::pow(d[26],2)*d[29]*std::pow(d[34],2) - d[2]*d[3]*std::pow(d[26],2)*d[29]*std::pow(d[35],2) - d[2]*d[3]*std::pow(d[26],2)*d[29]*std::pow(d[36],2) + d[2]*d[3]*std::pow(d[26],2)*d[29]*std::pow(d[37],2) - 2*d[2]*d[3]*d[26]*std::pow(d[28],2)*d[34]*d[37] - 2*d[2]*d[3]*d[26]*std::pow(d[28],2)*d[35]*d[36] - 2*d[2]*d[3]*d[26]*std::pow(d[29],2)*d[34]*d[37] - 2*d[2]*d[3]*d[26]*std::pow(d[29],2)*d[35]*d[36] + d[2]*d[3]*std::pow(d[27],2)*d[29]*std::pow(d[34],2) - d[2]*d[3]*std::pow(d[27],2)*d[29]*std::pow(d[35],2) - d[2]*d[3]*std::pow(d[27],2)*d[29]*std::pow(d[36],2) + d[2]*d[3]*std::pow(d[27],2)*d[29]*std::pow(d[37],2) - d[2]*d[3]*d[27]*std::pow(d[28],2)*std::pow(d[34],2) + d[2]*d[3]*d[27]*std::pow(d[28],2)*std::pow(d[35],2) - d[2]*d[3]*d[27]*std::pow(d[28],2)*std::pow(d[36],2) + d[2]*d[3]*d[27]*std::pow(d[28],2)*std::pow(d[37],2) - d[2]*d[3]*d[27]*std::pow(d[29],2)*std::pow(d[34],2) + d[2]*d[3]*d[27]*std::pow(d[29],2)*std::pow(d[35],2) - d[2]*d[3]*d[27]*std::pow(d[29],2)*std::pow(d[36],2) + d[2]*d[3]*d[27]*std::pow(d[29],2)*std::pow(d[37],2) + d[2]*d[21]*d[26]*d[38] + d[2]*d[21]*d[27]*d[39] - 2*d[2]*d[26]*d[34]*d[37] - 2*d[2]*d[26]*d[35]*d[36] - d[2]*d[27]*std::pow(d[34],2) + d[2]*d[27]*std::pow(d[35],2) - d[2]*d[27]*std::pow(d[36],2) + d[2]*d[27]*std::pow(d[37],2) + d[2]*d[29]*std::pow(d[34],2) - d[2]*d[29]*std::pow(d[35],2) - d[2]*d[29]*std::pow(d[36],2) + d[2]*d[29]*std::pow(d[37],2) - 2*d[26]*d[29]*d[34]*d[36] + 2*d[26]*d[29]*d[35]*d[37] + 2*d[27]*d[29]*d[34]*d[35] + 2*d[27]*d[29]*d[36]*d[37],  // NOLINT
        d[2]*std::pow(d[3],2)*d[5] - d[2]*d[3]*d[4] + d[2]*d[3]*d[8] - d[2]*d[7] + d[3]*d[11] - d[10], 0, d[2]*std::pow(d[3],2)*d[14] - d[2]*d[3]*d[13] + d[2]*d[3]*d[17] - d[2]*d[16] + d[3]*d[20] - d[19],  // NOLINT
        std::pow(d[2],2)*d[3]*d[6] + std::pow(d[2],2)*d[9] - d[2]*d[3]*d[5] - d[2]*d[8] + d[2]*d[12] - d[11], 0, std::pow(d[2],2)*d[3]*d[15] + std::pow(d[2],2)*d[18] - d[2]*d[3]*d[14] - d[2]*d[17] + d[2]*d[21] - d[20],  // NOLINT
        std::pow(d[2],2)*std::pow(d[3],2)*d[6] + std::pow(d[2],2)*d[3]*d[9] - d[2]*d[3]*d[4] + d[2]*d[3]*d[12] - d[2]*d[7] - d[10], 0, std::pow(d[2],2)*std::pow(d[3],2)*d[15] + std::pow(d[2],2)*d[3]*d[18] - d[2]*d[3]*d[13] + d[2]*d[3]*d[21] - d[2]*d[16] - d[19];  // NOLINT

        // Find the right eigenvectors of M2
        Eigen::JacobiSVD<Eigen::Matrix<double, 5, 3>> svd(M2, Eigen::ComputeFullV);

        // Extract coeffs
        d.head(2) = svd.matrixV().col(2).hnormalized();

        // Extract t2
        Eigen::Matrix<double, 5, 1> sols;
        sols << d(0), 0, d(1), d(2), d(3);
        Eigen::Matrix<double, 6, 1>  monoms;
        monoms << d(0) * d(2) * d(3), d(0) * d(2), d(0), d(2) * d(3), d(2), 1;
        sols(1) = -M.row(2) * monoms;

        return sols;
    }

    // TODO(marcusvaltonen): Move to helpers when necessary
    inline Eigen::Vector4d rot2quat(const Eigen::Matrix3d& R) {
        Eigen::Matrix4d K;
        K << R(0, 0)+R(1, 1)+R(2, 2), -R(1, 2)+R(2, 1), -R(2, 0)+R(0, 2), -R(0, 1)+R(1, 0),
            -R(1, 2)+R(2, 1), R(0, 0)-R(1, 1)-R(2, 2), R(1, 0)+R(0, 1), R(2, 0)+R(0, 2),
            -R(2, 0)+R(0, 2), R(1, 0)+R(0, 1), R(1, 1)-R(0, 0)-R(2, 2), R(2, 1)+R(1, 2),
            -R(0, 1)+R(1, 0), R(2, 0)+R(0, 2), R(2, 1)+R(1, 2), R(2, 2)-R(0, 0)-R(1, 1);

        K -= 3*Eigen::Matrix4d::Identity(4, 4);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(K, Eigen::ComputeFullV);

        // Extract coeffs
        Eigen::Vector4d v = svd.matrixV().col(3);
        return v;
    }

    // Function that utilizes the last equation of the DLT system to discard false solutions
    inline double get_algebraic_error_norot_frHfr(const Eigen::VectorXd& x) {
        // Compute algebraic error
        double error;
        double t2 = x(4)*x(4);
        double t3 = x(13)*x(13);
        double t4 = x(14)*x(14);
        double t5 = x(15)*x(15);
        double t6 = x(16)*x(16);
        double t7 = 1.0/x(3);
        double u1 = (x(1)*x(28) + x(24))*x(3);
        double u2 = t5 + t6;
        double u3 = t3 + t4;
        double w1 = x(28)*x(2) + x(25);
        double w2 = x(1)*x(26) + x(18);
        double w3 = x(1)*x(27) + x(21);

        error = (u3*w1*x(4) + t7*(x(26)*x(2) + x(19))*x(13) + t7*(x(27)*x(2) + x(22))*x(14) + w1)*x(16)
                + ((-u2 - u3)*x(4) - 1 - t2*(t5 + t6)*u3)*u1 - (u2*x(4) + 1)*(w2*x(13) + w3*x(14));

        return abs(error);
    }
}  // namespace HomLib::ValtonenOrnhagArxiv2020B
