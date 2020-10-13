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

#include "get_kukelova_cvpr_2015.hpp"
#include <Eigen/Geometry>
#include <vector>
#include "solver_kukelova_cvpr_2015.hpp"
#include "normalize2dpts.hpp"
#include "posedata.hpp"
#include "gj.hpp"

namespace HomLib::KukelovaCVPR2015 {
    inline Eigen::Matrix3d construct_homography_from_sols(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& tmp,
        const Eigen::MatrixXd& N);

    std::vector<HomLib::PoseData> get(const Eigen::MatrixXd &p1, const Eigen::MatrixXd &p2) {
        // This is a five point method
        const int nbr_pts = 5;

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

        // Compute distance to center for first points
        Eigen::VectorXd r21 = u1.colwise().squaredNorm();

        // Setup matrix for null space computation
        Eigen::MatrixXd M1(nbr_pts, 8);
        M1.setZero();

        for (int k = 0; k < nbr_pts; k++) {
            M1.row(k) << -r21(k) * u2(1, k), r21(k) * u2(0, k), -u2(1, k) * u1.col(k).homogeneous().transpose(),
                         u2(0, k) * u1.col(k).homogeneous().transpose();
        }

        // Find the null space
        // TODO(marcusvaltonen): This might be expensive - find out which is faster!
        // FullPivLU<MatrixXd> lu(M1);
        // MatrixXd N = lu.kernel();

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(M1, Eigen::ComputeFullV);
        Eigen::MatrixXd N = svd.matrixV().rightCols(3);

        // Create temporary input vector
        Eigen::MatrixXd tmp(4, nbr_pts);
        tmp << u1, u2;
        Eigen::VectorXd d(44);
        d << Eigen::Map<Eigen::VectorXd>(N.data(), 8*3),
             Eigen::Map<Eigen::VectorXd>(tmp.data(), 4*nbr_pts);

        // Create matrix M
        Eigen::MatrixXd M(7, 13);

        M << -d[24]*d[26], -d[25]*d[26], -std::pow(d[24],2)*d[26] - std::pow(d[25],2)*d[26], -d[26], d[0]*std::pow(d[24],2)*std::pow(d[26],2) + d[0]*std::pow(d[24],2)*std::pow(d[27],2) + d[0]*std::pow(d[25],2)*std::pow(d[26],2) + d[0]*std::pow(d[25],2)*std::pow(d[27],2) + d[2]*d[24]*std::pow(d[26],2) + d[2]*d[24]*std::pow(d[27],2) + d[3]*d[25]*std::pow(d[26],2) + d[3]*d[25]*std::pow(d[27],2) + d[4]*std::pow(d[26],2) + d[4]*std::pow(d[27],2), 0, d[0]*std::pow(d[24],2) + d[0]*std::pow(d[25],2) + d[2]*d[24] + d[3]*d[25] + d[4], d[8]*std::pow(d[24],2)*std::pow(d[26],2) + d[8]*std::pow(d[24],2)*std::pow(d[27],2) + d[8]*std::pow(d[25],2)*std::pow(d[26],2) + d[8]*std::pow(d[25],2)*std::pow(d[27],2) + d[10]*d[24]*std::pow(d[26],2) + d[10]*d[24]*std::pow(d[27],2) + d[11]*d[25]*std::pow(d[26],2) + d[11]*d[25]*std::pow(d[27],2) + d[12]*std::pow(d[26],2) + d[12]*std::pow(d[27],2), 0, d[16]*std::pow(d[24],2)*std::pow(d[26],2) + d[16]*std::pow(d[24],2)*std::pow(d[27],2) + d[16]*std::pow(d[25],2)*std::pow(d[26],2) + d[16]*std::pow(d[25],2)*std::pow(d[27],2) + d[18]*d[24]*std::pow(d[26],2) + d[18]*d[24]*std::pow(d[27],2) + d[19]*d[25]*std::pow(d[26],2) + d[19]*d[25]*std::pow(d[27],2) + d[20]*std::pow(d[26],2) + d[20]*std::pow(d[27],2), 0, d[8]*std::pow(d[24],2) + d[8]*std::pow(d[25],2) + d[10]*d[24] + d[11]*d[25] + d[12], d[16]*std::pow(d[24],2) + d[16]*std::pow(d[25],2) + d[18]*d[24] + d[19]*d[25] + d[20],  // NOLINT
        -d[28]*d[30], -d[29]*d[30], -std::pow(d[28],2)*d[30] - std::pow(d[29],2)*d[30], -d[30], d[0]*std::pow(d[28],2)*std::pow(d[30],2) + d[0]*std::pow(d[28],2)*std::pow(d[31],2) + d[0]*std::pow(d[29],2)*std::pow(d[30],2) + d[0]*std::pow(d[29],2)*std::pow(d[31],2) + d[2]*d[28]*std::pow(d[30],2) + d[2]*d[28]*std::pow(d[31],2) + d[3]*d[29]*std::pow(d[30],2) + d[3]*d[29]*std::pow(d[31],2) + d[4]*std::pow(d[30],2) + d[4]*std::pow(d[31],2), 0, d[0]*std::pow(d[28],2) + d[0]*std::pow(d[29],2) + d[2]*d[28] + d[3]*d[29] + d[4], d[8]*std::pow(d[28],2)*std::pow(d[30],2) + d[8]*std::pow(d[28],2)*std::pow(d[31],2) + d[8]*std::pow(d[29],2)*std::pow(d[30],2) + d[8]*std::pow(d[29],2)*std::pow(d[31],2) + d[10]*d[28]*std::pow(d[30],2) + d[10]*d[28]*std::pow(d[31],2) + d[11]*d[29]*std::pow(d[30],2) + d[11]*d[29]*std::pow(d[31],2) + d[12]*std::pow(d[30],2) + d[12]*std::pow(d[31],2), 0, d[16]*std::pow(d[28],2)*std::pow(d[30],2) + d[16]*std::pow(d[28],2)*std::pow(d[31],2) + d[16]*std::pow(d[29],2)*std::pow(d[30],2) + d[16]*std::pow(d[29],2)*std::pow(d[31],2) + d[18]*d[28]*std::pow(d[30],2) + d[18]*d[28]*std::pow(d[31],2) + d[19]*d[29]*std::pow(d[30],2) + d[19]*d[29]*std::pow(d[31],2) + d[20]*std::pow(d[30],2) + d[20]*std::pow(d[31],2), 0, d[8]*std::pow(d[28],2) + d[8]*std::pow(d[29],2) + d[10]*d[28] + d[11]*d[29] + d[12], d[16]*std::pow(d[28],2) + d[16]*std::pow(d[29],2) + d[18]*d[28] + d[19]*d[29] + d[20],  // NOLINT
        -d[32]*d[34], -d[33]*d[34], -std::pow(d[32],2)*d[34] - std::pow(d[33],2)*d[34], -d[34], d[0]*std::pow(d[32],2)*std::pow(d[34],2) + d[0]*std::pow(d[32],2)*std::pow(d[35],2) + d[0]*std::pow(d[33],2)*std::pow(d[34],2) + d[0]*std::pow(d[33],2)*std::pow(d[35],2) + d[2]*d[32]*std::pow(d[34],2) + d[2]*d[32]*std::pow(d[35],2) + d[3]*d[33]*std::pow(d[34],2) + d[3]*d[33]*std::pow(d[35],2) + d[4]*std::pow(d[34],2) + d[4]*std::pow(d[35],2), 0, d[0]*std::pow(d[32],2) + d[0]*std::pow(d[33],2) + d[2]*d[32] + d[3]*d[33] + d[4], d[8]*std::pow(d[32],2)*std::pow(d[34],2) + d[8]*std::pow(d[32],2)*std::pow(d[35],2) + d[8]*std::pow(d[33],2)*std::pow(d[34],2) + d[8]*std::pow(d[33],2)*std::pow(d[35],2) + d[10]*d[32]*std::pow(d[34],2) + d[10]*d[32]*std::pow(d[35],2) + d[11]*d[33]*std::pow(d[34],2) + d[11]*d[33]*std::pow(d[35],2) + d[12]*std::pow(d[34],2) + d[12]*std::pow(d[35],2), 0, d[16]*std::pow(d[32],2)*std::pow(d[34],2) + d[16]*std::pow(d[32],2)*std::pow(d[35],2) + d[16]*std::pow(d[33],2)*std::pow(d[34],2) + d[16]*std::pow(d[33],2)*std::pow(d[35],2) + d[18]*d[32]*std::pow(d[34],2) + d[18]*d[32]*std::pow(d[35],2) + d[19]*d[33]*std::pow(d[34],2) + d[19]*d[33]*std::pow(d[35],2) + d[20]*std::pow(d[34],2) + d[20]*std::pow(d[35],2), 0, d[8]*std::pow(d[32],2) + d[8]*std::pow(d[33],2) + d[10]*d[32] + d[11]*d[33] + d[12], d[16]*std::pow(d[32],2) + d[16]*std::pow(d[33],2) + d[18]*d[32] + d[19]*d[33] + d[20],  // NOLINT
        -d[36]*d[38], -d[37]*d[38], -std::pow(d[36],2)*d[38] - std::pow(d[37],2)*d[38], -d[38], d[0]*std::pow(d[36],2)*std::pow(d[38],2) + d[0]*std::pow(d[36],2)*std::pow(d[39],2) + d[0]*std::pow(d[37],2)*std::pow(d[38],2) + d[0]*std::pow(d[37],2)*std::pow(d[39],2) + d[2]*d[36]*std::pow(d[38],2) + d[2]*d[36]*std::pow(d[39],2) + d[3]*d[37]*std::pow(d[38],2) + d[3]*d[37]*std::pow(d[39],2) + d[4]*std::pow(d[38],2) + d[4]*std::pow(d[39],2), 0, d[0]*std::pow(d[36],2) + d[0]*std::pow(d[37],2) + d[2]*d[36] + d[3]*d[37] + d[4], d[8]*std::pow(d[36],2)*std::pow(d[38],2) + d[8]*std::pow(d[36],2)*std::pow(d[39],2) + d[8]*std::pow(d[37],2)*std::pow(d[38],2) + d[8]*std::pow(d[37],2)*std::pow(d[39],2) + d[10]*d[36]*std::pow(d[38],2) + d[10]*d[36]*std::pow(d[39],2) + d[11]*d[37]*std::pow(d[38],2) + d[11]*d[37]*std::pow(d[39],2) + d[12]*std::pow(d[38],2) + d[12]*std::pow(d[39],2), 0, d[16]*std::pow(d[36],2)*std::pow(d[38],2) + d[16]*std::pow(d[36],2)*std::pow(d[39],2) + d[16]*std::pow(d[37],2)*std::pow(d[38],2) + d[16]*std::pow(d[37],2)*std::pow(d[39],2) + d[18]*d[36]*std::pow(d[38],2) + d[18]*d[36]*std::pow(d[39],2) + d[19]*d[37]*std::pow(d[38],2) + d[19]*d[37]*std::pow(d[39],2) + d[20]*std::pow(d[38],2) + d[20]*std::pow(d[39],2), 0, d[8]*std::pow(d[36],2) + d[8]*std::pow(d[37],2) + d[10]*d[36] + d[11]*d[37] + d[12], d[16]*std::pow(d[36],2) + d[16]*std::pow(d[37],2) + d[18]*d[36] + d[19]*d[37] + d[20],  // NOLINT
        -d[40]*d[42], -d[41]*d[42], -std::pow(d[40],2)*d[42] - std::pow(d[41],2)*d[42], -d[42], d[0]*std::pow(d[40],2)*std::pow(d[42],2) + d[0]*std::pow(d[40],2)*std::pow(d[43],2) + d[0]*std::pow(d[41],2)*std::pow(d[42],2) + d[0]*std::pow(d[41],2)*std::pow(d[43],2) + d[2]*d[40]*std::pow(d[42],2) + d[2]*d[40]*std::pow(d[43],2) + d[3]*d[41]*std::pow(d[42],2) + d[3]*d[41]*std::pow(d[43],2) + d[4]*std::pow(d[42],2) + d[4]*std::pow(d[43],2), 0, d[0]*std::pow(d[40],2) + d[0]*std::pow(d[41],2) + d[2]*d[40] + d[3]*d[41] + d[4], d[8]*std::pow(d[40],2)*std::pow(d[42],2) + d[8]*std::pow(d[40],2)*std::pow(d[43],2) + d[8]*std::pow(d[41],2)*std::pow(d[42],2) + d[8]*std::pow(d[41],2)*std::pow(d[43],2) + d[10]*d[40]*std::pow(d[42],2) + d[10]*d[40]*std::pow(d[43],2) + d[11]*d[41]*std::pow(d[42],2) + d[11]*d[41]*std::pow(d[43],2) + d[12]*std::pow(d[42],2) + d[12]*std::pow(d[43],2), 0, d[16]*std::pow(d[40],2)*std::pow(d[42],2) + d[16]*std::pow(d[40],2)*std::pow(d[43],2) + d[16]*std::pow(d[41],2)*std::pow(d[42],2) + d[16]*std::pow(d[41],2)*std::pow(d[43],2) + d[18]*d[40]*std::pow(d[42],2) + d[18]*d[40]*std::pow(d[43],2) + d[19]*d[41]*std::pow(d[42],2) + d[19]*d[41]*std::pow(d[43],2) + d[20]*std::pow(d[42],2) + d[20]*std::pow(d[43],2), 0, d[8]*std::pow(d[40],2) + d[8]*std::pow(d[41],2) + d[10]*d[40] + d[11]*d[41] + d[12], d[16]*std::pow(d[40],2) + d[16]*std::pow(d[41],2) + d[18]*d[40] + d[19]*d[41] + d[20],  // NOLINT
        0, 0, 0, 0, 0, -d[4], d[0], 0, -d[12], 0, -d[20], d[8], d[16],
        0, 0, 0, 0, 0, -d[7], d[1], 0, -d[15], 0, -d[23], d[9], d[17];

        // Gauss-Jordan
        HomLib::gj(&M);

        // Wrap input data to expected format
        Eigen::VectorXd input(66);
        input << Eigen::Map<Eigen::VectorXd>(N.data(), 8*3),
                 Eigen::Map<Eigen::VectorXd>(M.rightCols(6).data(), 6*7);

        // Extract solution
        Eigen::MatrixXcd sols = solver_kukelova_cvpr_2015(input);

        // Pre-processing: Remove complex-valued solutions
        double thresh = 1e-5;
        Eigen::ArrayXd real_sols(3);
        real_sols = sols.imag().cwiseAbs().colwise().sum();
        int nbr_real_sols = (real_sols <= thresh).count();

        // Create putative solutions
        Eigen::Matrix3d Htmp;
        std::complex<double> lam1;
        std::complex<double> lam2;
        std::vector<HomLib::PoseData> posedata(nbr_real_sols);
        Eigen::ArrayXd xx(3);

        int cnt = 0;
        for (int i = 0; i < real_sols.size(); i++) {
            if (real_sols(i) <= thresh) {
                // Get parameters.
                xx = sols.col(i).real();

                // Construct putative fundamental matrix
                Htmp = S.inverse() * HomLib::KukelovaCVPR2015::construct_homography_from_sols(xx, d, N) * S;

                // Package output
                posedata[cnt].homography = Htmp;
                posedata[cnt].distortion_parameter = xx(0) * std::pow(scale, 2);
                posedata[cnt].distortion_parameter2 = xx(1) * std::pow(scale, 2);
                cnt++;
            }
        }

        return posedata;
    }

    // Function that utilizes the last equation of the DLT system to discard false solutions
    inline Eigen::Matrix3d construct_homography_from_sols(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& tmp,
        const Eigen::MatrixXd& N
    ) {
        Eigen::VectorXd d(51);
        d << 0, 0, 0, xx(0), xx(1), 0, xx(2), tmp;

        Eigen::Matrix<double, 5, 5> M;
        M << -d[31]*d[33], -d[32]*d[33], -d[3]*std::pow(d[31],2)*d[33] - d[3]*std::pow(d[32],2)*d[33] - d[33], d[4]*d[7]*std::pow(d[31],2)*std::pow(d[33],2) + d[4]*d[7]*std::pow(d[31],2)*std::pow(d[34],2) + d[4]*d[7]*std::pow(d[32],2)*std::pow(d[33],2) + d[4]*d[7]*std::pow(d[32],2)*std::pow(d[34],2) + d[4]*d[9]*d[31]*std::pow(d[33],2) + d[4]*d[9]*d[31]*std::pow(d[34],2) + d[4]*d[10]*d[32]*std::pow(d[33],2) + d[4]*d[10]*d[32]*std::pow(d[34],2) + d[4]*d[11]*std::pow(d[33],2) + d[4]*d[11]*std::pow(d[34],2) + d[7]*std::pow(d[31],2) + d[7]*std::pow(d[32],2) + d[9]*d[31] + d[10]*d[32] + d[11], d[4]*d[6]*d[15]*std::pow(d[31],2)*std::pow(d[33],2) + d[4]*d[6]*d[15]*std::pow(d[31],2)*std::pow(d[34],2) + d[4]*d[6]*d[15]*std::pow(d[32],2)*std::pow(d[33],2) + d[4]*d[6]*d[15]*std::pow(d[32],2)*std::pow(d[34],2) + d[4]*d[6]*d[17]*d[31]*std::pow(d[33],2) + d[4]*d[6]*d[17]*d[31]*std::pow(d[34],2) + d[4]*d[6]*d[18]*d[32]*std::pow(d[33],2) + d[4]*d[6]*d[18]*d[32]*std::pow(d[34],2) + d[4]*d[6]*d[19]*std::pow(d[33],2) + d[4]*d[6]*d[19]*std::pow(d[34],2) + d[4]*d[23]*std::pow(d[31],2)*std::pow(d[33],2) + d[4]*d[23]*std::pow(d[31],2)*std::pow(d[34],2) + d[4]*d[23]*std::pow(d[32],2)*std::pow(d[33],2) + d[4]*d[23]*std::pow(d[32],2)*std::pow(d[34],2) + d[4]*d[25]*d[31]*std::pow(d[33],2) + d[4]*d[25]*d[31]*std::pow(d[34],2) + d[4]*d[26]*d[32]*std::pow(d[33],2) + d[4]*d[26]*d[32]*std::pow(d[34],2) + d[4]*d[27]*std::pow(d[33],2) + d[4]*d[27]*std::pow(d[34],2) + d[6]*d[15]*std::pow(d[31],2) + d[6]*d[15]*std::pow(d[32],2) + d[6]*d[17]*d[31] + d[6]*d[18]*d[32] + d[6]*d[19] + d[23]*std::pow(d[31],2) + d[23]*std::pow(d[32],2) + d[25]*d[31] + d[26]*d[32] + d[27],  // NOLINT
        -d[35]*d[37], -d[36]*d[37], -d[3]*std::pow(d[35],2)*d[37] - d[3]*std::pow(d[36],2)*d[37] - d[37], d[4]*d[7]*std::pow(d[35],2)*std::pow(d[37],2) + d[4]*d[7]*std::pow(d[35],2)*std::pow(d[38],2) + d[4]*d[7]*std::pow(d[36],2)*std::pow(d[37],2) + d[4]*d[7]*std::pow(d[36],2)*std::pow(d[38],2) + d[4]*d[9]*d[35]*std::pow(d[37],2) + d[4]*d[9]*d[35]*std::pow(d[38],2) + d[4]*d[10]*d[36]*std::pow(d[37],2) + d[4]*d[10]*d[36]*std::pow(d[38],2) + d[4]*d[11]*std::pow(d[37],2) + d[4]*d[11]*std::pow(d[38],2) + d[7]*std::pow(d[35],2) + d[7]*std::pow(d[36],2) + d[9]*d[35] + d[10]*d[36] + d[11], d[4]*d[6]*d[15]*std::pow(d[35],2)*std::pow(d[37],2) + d[4]*d[6]*d[15]*std::pow(d[35],2)*std::pow(d[38],2) + d[4]*d[6]*d[15]*std::pow(d[36],2)*std::pow(d[37],2) + d[4]*d[6]*d[15]*std::pow(d[36],2)*std::pow(d[38],2) + d[4]*d[6]*d[17]*d[35]*std::pow(d[37],2) + d[4]*d[6]*d[17]*d[35]*std::pow(d[38],2) + d[4]*d[6]*d[18]*d[36]*std::pow(d[37],2) + d[4]*d[6]*d[18]*d[36]*std::pow(d[38],2) + d[4]*d[6]*d[19]*std::pow(d[37],2) + d[4]*d[6]*d[19]*std::pow(d[38],2) + d[4]*d[23]*std::pow(d[35],2)*std::pow(d[37],2) + d[4]*d[23]*std::pow(d[35],2)*std::pow(d[38],2) + d[4]*d[23]*std::pow(d[36],2)*std::pow(d[37],2) + d[4]*d[23]*std::pow(d[36],2)*std::pow(d[38],2) + d[4]*d[25]*d[35]*std::pow(d[37],2) + d[4]*d[25]*d[35]*std::pow(d[38],2) + d[4]*d[26]*d[36]*std::pow(d[37],2) + d[4]*d[26]*d[36]*std::pow(d[38],2) + d[4]*d[27]*std::pow(d[37],2) + d[4]*d[27]*std::pow(d[38],2) + d[6]*d[15]*std::pow(d[35],2) + d[6]*d[15]*std::pow(d[36],2) + d[6]*d[17]*d[35] + d[6]*d[18]*d[36] + d[6]*d[19] + d[23]*std::pow(d[35],2) + d[23]*std::pow(d[36],2) + d[25]*d[35] + d[26]*d[36] + d[27],  // NOLINT
        -d[39]*d[41], -d[40]*d[41], -d[3]*std::pow(d[39],2)*d[41] - d[3]*std::pow(d[40],2)*d[41] - d[41], d[4]*d[7]*std::pow(d[39],2)*std::pow(d[41],2) + d[4]*d[7]*std::pow(d[39],2)*std::pow(d[42],2) + d[4]*d[7]*std::pow(d[40],2)*std::pow(d[41],2) + d[4]*d[7]*std::pow(d[40],2)*std::pow(d[42],2) + d[4]*d[9]*d[39]*std::pow(d[41],2) + d[4]*d[9]*d[39]*std::pow(d[42],2) + d[4]*d[10]*d[40]*std::pow(d[41],2) + d[4]*d[10]*d[40]*std::pow(d[42],2) + d[4]*d[11]*std::pow(d[41],2) + d[4]*d[11]*std::pow(d[42],2) + d[7]*std::pow(d[39],2) + d[7]*std::pow(d[40],2) + d[9]*d[39] + d[10]*d[40] + d[11], d[4]*d[6]*d[15]*std::pow(d[39],2)*std::pow(d[41],2) + d[4]*d[6]*d[15]*std::pow(d[39],2)*std::pow(d[42],2) + d[4]*d[6]*d[15]*std::pow(d[40],2)*std::pow(d[41],2) + d[4]*d[6]*d[15]*std::pow(d[40],2)*std::pow(d[42],2) + d[4]*d[6]*d[17]*d[39]*std::pow(d[41],2) + d[4]*d[6]*d[17]*d[39]*std::pow(d[42],2) + d[4]*d[6]*d[18]*d[40]*std::pow(d[41],2) + d[4]*d[6]*d[18]*d[40]*std::pow(d[42],2) + d[4]*d[6]*d[19]*std::pow(d[41],2) + d[4]*d[6]*d[19]*std::pow(d[42],2) + d[4]*d[23]*std::pow(d[39],2)*std::pow(d[41],2) + d[4]*d[23]*std::pow(d[39],2)*std::pow(d[42],2) + d[4]*d[23]*std::pow(d[40],2)*std::pow(d[41],2) + d[4]*d[23]*std::pow(d[40],2)*std::pow(d[42],2) + d[4]*d[25]*d[39]*std::pow(d[41],2) + d[4]*d[25]*d[39]*std::pow(d[42],2) + d[4]*d[26]*d[40]*std::pow(d[41],2) + d[4]*d[26]*d[40]*std::pow(d[42],2) + d[4]*d[27]*std::pow(d[41],2) + d[4]*d[27]*std::pow(d[42],2) + d[6]*d[15]*std::pow(d[39],2) + d[6]*d[15]*std::pow(d[40],2) + d[6]*d[17]*d[39] + d[6]*d[18]*d[40] + d[6]*d[19] + d[23]*std::pow(d[39],2) + d[23]*std::pow(d[40],2) + d[25]*d[39] + d[26]*d[40] + d[27],  // NOLINT
        -d[43]*d[45], -d[44]*d[45], -d[3]*std::pow(d[43],2)*d[45] - d[3]*std::pow(d[44],2)*d[45] - d[45], d[4]*d[7]*std::pow(d[43],2)*std::pow(d[45],2) + d[4]*d[7]*std::pow(d[43],2)*std::pow(d[46],2) + d[4]*d[7]*std::pow(d[44],2)*std::pow(d[45],2) + d[4]*d[7]*std::pow(d[44],2)*std::pow(d[46],2) + d[4]*d[9]*d[43]*std::pow(d[45],2) + d[4]*d[9]*d[43]*std::pow(d[46],2) + d[4]*d[10]*d[44]*std::pow(d[45],2) + d[4]*d[10]*d[44]*std::pow(d[46],2) + d[4]*d[11]*std::pow(d[45],2) + d[4]*d[11]*std::pow(d[46],2) + d[7]*std::pow(d[43],2) + d[7]*std::pow(d[44],2) + d[9]*d[43] + d[10]*d[44] + d[11], d[4]*d[6]*d[15]*std::pow(d[43],2)*std::pow(d[45],2) + d[4]*d[6]*d[15]*std::pow(d[43],2)*std::pow(d[46],2) + d[4]*d[6]*d[15]*std::pow(d[44],2)*std::pow(d[45],2) + d[4]*d[6]*d[15]*std::pow(d[44],2)*std::pow(d[46],2) + d[4]*d[6]*d[17]*d[43]*std::pow(d[45],2) + d[4]*d[6]*d[17]*d[43]*std::pow(d[46],2) + d[4]*d[6]*d[18]*d[44]*std::pow(d[45],2) + d[4]*d[6]*d[18]*d[44]*std::pow(d[46],2) + d[4]*d[6]*d[19]*std::pow(d[45],2) + d[4]*d[6]*d[19]*std::pow(d[46],2) + d[4]*d[23]*std::pow(d[43],2)*std::pow(d[45],2) + d[4]*d[23]*std::pow(d[43],2)*std::pow(d[46],2) + d[4]*d[23]*std::pow(d[44],2)*std::pow(d[45],2) + d[4]*d[23]*std::pow(d[44],2)*std::pow(d[46],2) + d[4]*d[25]*d[43]*std::pow(d[45],2) + d[4]*d[25]*d[43]*std::pow(d[46],2) + d[4]*d[26]*d[44]*std::pow(d[45],2) + d[4]*d[26]*d[44]*std::pow(d[46],2) + d[4]*d[27]*std::pow(d[45],2) + d[4]*d[27]*std::pow(d[46],2) + d[6]*d[15]*std::pow(d[43],2) + d[6]*d[15]*std::pow(d[44],2) + d[6]*d[17]*d[43] + d[6]*d[18]*d[44] + d[6]*d[19] + d[23]*std::pow(d[43],2) + d[23]*std::pow(d[44],2) + d[25]*d[43] + d[26]*d[44] + d[27],  // NOLINT
        -d[47]*d[49], -d[48]*d[49], -d[3]*std::pow(d[47],2)*d[49] - d[3]*std::pow(d[48],2)*d[49] - d[49], d[4]*d[7]*std::pow(d[47],2)*std::pow(d[49],2) + d[4]*d[7]*std::pow(d[47],2)*std::pow(d[50],2) + d[4]*d[7]*std::pow(d[48],2)*std::pow(d[49],2) + d[4]*d[7]*std::pow(d[48],2)*std::pow(d[50],2) + d[4]*d[9]*d[47]*std::pow(d[49],2) + d[4]*d[9]*d[47]*std::pow(d[50],2) + d[4]*d[10]*d[48]*std::pow(d[49],2) + d[4]*d[10]*d[48]*std::pow(d[50],2) + d[4]*d[11]*std::pow(d[49],2) + d[4]*d[11]*std::pow(d[50],2) + d[7]*std::pow(d[47],2) + d[7]*std::pow(d[48],2) + d[9]*d[47] + d[10]*d[48] + d[11], d[4]*d[6]*d[15]*std::pow(d[47],2)*std::pow(d[49],2) + d[4]*d[6]*d[15]*std::pow(d[47],2)*std::pow(d[50],2) + d[4]*d[6]*d[15]*std::pow(d[48],2)*std::pow(d[49],2) + d[4]*d[6]*d[15]*std::pow(d[48],2)*std::pow(d[50],2) + d[4]*d[6]*d[17]*d[47]*std::pow(d[49],2) + d[4]*d[6]*d[17]*d[47]*std::pow(d[50],2) + d[4]*d[6]*d[18]*d[48]*std::pow(d[49],2) + d[4]*d[6]*d[18]*d[48]*std::pow(d[50],2) + d[4]*d[6]*d[19]*std::pow(d[49],2) + d[4]*d[6]*d[19]*std::pow(d[50],2) + d[4]*d[23]*std::pow(d[47],2)*std::pow(d[49],2) + d[4]*d[23]*std::pow(d[47],2)*std::pow(d[50],2) + d[4]*d[23]*std::pow(d[48],2)*std::pow(d[49],2) + d[4]*d[23]*std::pow(d[48],2)*std::pow(d[50],2) + d[4]*d[25]*d[47]*std::pow(d[49],2) + d[4]*d[25]*d[47]*std::pow(d[50],2) + d[4]*d[26]*d[48]*std::pow(d[49],2) + d[4]*d[26]*d[48]*std::pow(d[50],2) + d[4]*d[27]*std::pow(d[49],2) + d[4]*d[27]*std::pow(d[50],2) + d[6]*d[15]*std::pow(d[47],2) + d[6]*d[15]*std::pow(d[48],2) + d[6]*d[17]*d[47] + d[6]*d[18]*d[48] + d[6]*d[19] + d[23]*std::pow(d[47],2) + d[23]*std::pow(d[48],2) + d[25]*d[47] + d[26]*d[48] + d[27];  // NOLINT

        // Perform SVD and extract coeffs
        Eigen::JacobiSVD<Eigen::Matrix<double, 5, 5>> svd(M, Eigen::ComputeFullV);
        Eigen::Matrix<double, 4, 1>  v = svd.matrixV().col(4).hnormalized();

        // Construct homography
        Eigen::Matrix<double, 8, 1> v1 = N.col(0) * v(3) + N.col(1) * xx(2) + N.col(2);
        Eigen::Matrix3d H;
        H.row(0) = v1.segment<3>(2);
        H.row(1) = v1.segment<3>(5);
        H.row(2) = v.head<3>();

        return H;
    }
}  // namespace HomLib::KukelovaCVPR2015
