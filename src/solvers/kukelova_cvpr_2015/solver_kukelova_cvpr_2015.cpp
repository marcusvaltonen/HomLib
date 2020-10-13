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

#include "solver_kukelova_cvpr_2015.hpp"
#include <Eigen/Dense>

namespace HomLib::KukelovaCVPR2015 {
    Eigen::MatrixXcd solver_kukelova_cvpr_2015(const Eigen::VectorXd& data) {
        // Compute coefficients
        const double* d = data.data();
        Eigen::VectorXd coeffs(22);
        coeffs[0] = -d[34];
        coeffs[1] = -d[27];
        coeffs[2] = -d[48];
        coeffs[3] = -d[41];
        coeffs[4] = d[33] - d[55];
        coeffs[5] = d[26];
        coeffs[6] = d[47] - d[62];
        coeffs[7] = d[40];
        coeffs[8] = d[54];
        coeffs[9] = d[61];
        coeffs[10] = -d[37];
        coeffs[11] = -d[51];
        coeffs[12] = d[35];
        coeffs[13] = d[28] - d[58];
        coeffs[14] = d[49];
        coeffs[15] = d[42] - d[65];
        coeffs[16] = d[56];
        coeffs[17] = d[63];
        coeffs[18] = d[36] - d[58];
        coeffs[19] = d[50] - d[65];
        coeffs[20] = d[57];
        coeffs[21] = d[64];

        // Setup elimination template
        static const int coeffs0_ind[] = { 10,10,0,10,1,10,10,12,2,0,10,11,13,3,11,1,10,10,18,11,4,12,0,10,18,5,13,1,10,18,15,3,11,11,19,16,6,14,4,12,2,18,0,10,11,19,7,15,5,13,3,11,1,18,20,10,19,5,13,20,7,15,5,20,13,21,17,6,14,19,2,11,7,15,3,19,21,11,9,17,8,16,6,14,20,4,18,19,12,21,8,16,4,12,18,20 };  // NOLINT
        static const int coeffs1_ind[] = { 9,21,17,9,17,21,6,19,14,7,21,15,9,17,8,20,21,16,8,16,20 };  // NOLINT
        static const int C0_ind[] = { 0,12,17,30,33,34,47,48,49,51,55,62,64,65,66,67,68,75,76,79,81,82,85,90,94,97,98,101,102,111,112,115,116,123,124,128,129,130,131,132,133,135,136,137,138,142,145,146,147,148,149,150,152,155,156,157,159,165,166,175,181,182,184,187,189,191,192,195,196,199,200,201,211,212,216,219,220,221,225,226,227,228,229,230,231,232,233,234,237,238,241,242,245,246,250,254 };  // NOLINT
        static const int C1_ind[] = { 8,9,13,19,20,23,24,25,29,40,43,45,53,54,56,57,58,61,69,70,74 };  // NOLINT

        Eigen::Matrix<double, 16, 16> C0; C0.setZero();
        Eigen::Matrix<double, 16, 5> C1; C1.setZero();
        for (int i = 0; i < 96; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
        for (int i = 0; i < 21; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); }

        Eigen::Matrix<double, 16, 5> C12 = C0.partialPivLu().solve(C1);

        // Setup action matrix
        Eigen::Matrix<double, 9, 5> RR;
        RR << -C12.bottomRows(4), Eigen::Matrix<double, 5, 5>::Identity(5, 5);

        static const int AM_ind[] = { 5,0,1,2,3 };  // NOLINT
        Eigen::Matrix<double, 5, 5> AM;
        for (int i = 0; i < 5; i++) {
            AM.row(i) = RR.row(AM_ind[i]);
        }

        Eigen::Matrix<std::complex<double>, 3, 5> sols;
        sols.setZero();

        // Solve eigenvalue problem
        Eigen::EigenSolver<Eigen::Matrix<double, 5, 5> > es(AM);
        Eigen::ArrayXcd D = es.eigenvalues();
        Eigen::ArrayXXcd V = es.eigenvectors();
        V = (V / V.row(0).array().replicate(5, 1)).eval();

        sols.row(0) = D.transpose().array();
        sols.row(1) = V.row(2).array();
        sols.row(2) = V.row(3).array();

        return sols;
    }
}  // namespace HomLib::KukelovaCVPR2015
