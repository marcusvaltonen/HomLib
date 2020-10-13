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

#include <Eigen/Dense>

#include "solver_valtonenornhag_arxiv_2020a_fHf.hpp"

namespace HomLib::ValtonenOrnhagArxiv2020A {
    Eigen::MatrixXcd solver_fHf(const Eigen::VectorXd& data) {
        // Compute coefficients
        const double* d = data.data();
        Eigen::VectorXd coeffs(45);
        coeffs[0] = 1;
        coeffs[1] = -1;
        coeffs[2] = d[20]*d[28];
        coeffs[3] = -d[18]*d[28];
        coeffs[4] = -d[19]*d[29];
        coeffs[5] = d[19]*d[28];
        coeffs[6] = d[2]*d[20]*d[22] + d[3]*d[20]*d[25] + d[0]*d[14]*d[28] + d[1]*d[17]*d[28];
        coeffs[7] = -d[2]*d[18]*d[22] - d[3]*d[18]*d[25] - d[0]*d[12]*d[28] - d[1]*d[15]*d[28];
        coeffs[8] = -d[2]*d[19]*d[23] - d[3]*d[19]*d[26] - d[0]*d[13]*d[29] - d[1]*d[16]*d[29];
        coeffs[9] = d[2]*d[19]*d[22] + d[3]*d[19]*d[25] + d[0]*d[13]*d[28] + d[1]*d[16]*d[28];
        coeffs[10] = d[0]*d[2]*d[14]*d[22] + d[1]*d[2]*d[17]*d[22] + d[0]*d[3]*d[14]*d[25] + d[1]*d[3]*d[17]*d[25];
        coeffs[11] = -d[0]*d[2]*d[12]*d[22] - d[1]*d[2]*d[15]*d[22] - d[0]*d[3]*d[12]*d[25] - d[1]*d[3]*d[15]*d[25];
        coeffs[12] = -d[0]*d[2]*d[13]*d[23] - d[1]*d[2]*d[16]*d[23] - d[0]*d[3]*d[13]*d[26] - d[1]*d[3]*d[16]*d[26];
        coeffs[13] = d[0]*d[2]*d[13]*d[22] + d[1]*d[2]*d[16]*d[22] + d[0]*d[3]*d[13]*d[25] + d[1]*d[3]*d[16]*d[25];
        coeffs[14] = -d[20]*d[28];
        coeffs[15] = -d[19]*d[28];
        coeffs[16] = d[19]*d[27];
        coeffs[17] = -d[2]*d[20]*d[22] - d[3]*d[20]*d[25] - d[0]*d[14]*d[28] - d[1]*d[17]*d[28];
        coeffs[18] = -d[2]*d[19]*d[22] - d[3]*d[19]*d[25] - d[0]*d[13]*d[28] - d[1]*d[16]*d[28];
        coeffs[19] = d[2]*d[19]*d[21] + d[3]*d[19]*d[24] + d[0]*d[13]*d[27] + d[1]*d[16]*d[27];
        coeffs[20] = -d[0]*d[2]*d[14]*d[22] - d[1]*d[2]*d[17]*d[22] - d[0]*d[3]*d[14]*d[25] - d[1]*d[3]*d[17]*d[25];
        coeffs[21] = -d[0]*d[2]*d[13]*d[22] - d[1]*d[2]*d[16]*d[22] - d[0]*d[3]*d[13]*d[25] - d[1]*d[3]*d[16]*d[25];
        coeffs[22] = d[0]*d[2]*d[13]*d[21] + d[1]*d[2]*d[16]*d[21] + d[0]*d[3]*d[13]*d[24] + d[1]*d[3]*d[16]*d[24];
        coeffs[23] = d[6]*d[20]*d[22] + d[7]*d[20]*d[25] + d[4]*d[14]*d[28] + d[5]*d[17]*d[28];
        coeffs[24] = -d[6]*d[18]*d[22] - d[7]*d[18]*d[25] - d[4]*d[12]*d[28] - d[5]*d[15]*d[28];
        coeffs[25] = -d[6]*d[19]*d[23] - d[7]*d[19]*d[26] - d[4]*d[13]*d[29] - d[5]*d[16]*d[29];
        coeffs[26] = d[6]*d[19]*d[22] + d[7]*d[19]*d[25] + d[4]*d[13]*d[28] + d[5]*d[16]*d[28];
        coeffs[27] = d[4]*d[6]*d[14]*d[22] + d[5]*d[6]*d[17]*d[22] + d[4]*d[7]*d[14]*d[25] + d[5]*d[7]*d[17]*d[25];
        coeffs[28] = -d[4]*d[6]*d[12]*d[22] - d[5]*d[6]*d[15]*d[22] - d[4]*d[7]*d[12]*d[25] - d[5]*d[7]*d[15]*d[25];
        coeffs[29] = -d[4]*d[6]*d[13]*d[23] - d[5]*d[6]*d[16]*d[23] - d[4]*d[7]*d[13]*d[26] - d[5]*d[7]*d[16]*d[26];
        coeffs[30] = d[4]*d[6]*d[13]*d[22] + d[5]*d[6]*d[16]*d[22] + d[4]*d[7]*d[13]*d[25] + d[5]*d[7]*d[16]*d[25];
        coeffs[31] = -d[6]*d[20]*d[22] - d[7]*d[20]*d[25] - d[4]*d[14]*d[28] - d[5]*d[17]*d[28];
        coeffs[32] = -d[6]*d[19]*d[22] - d[7]*d[19]*d[25] - d[4]*d[13]*d[28] - d[5]*d[16]*d[28];
        coeffs[33] = d[6]*d[19]*d[21] + d[7]*d[19]*d[24] + d[4]*d[13]*d[27] + d[5]*d[16]*d[27];
        coeffs[34] = -d[4]*d[6]*d[14]*d[22] - d[5]*d[6]*d[17]*d[22] - d[4]*d[7]*d[14]*d[25] - d[5]*d[7]*d[17]*d[25];
        coeffs[35] = -d[4]*d[6]*d[13]*d[22] - d[5]*d[6]*d[16]*d[22] - d[4]*d[7]*d[13]*d[25] - d[5]*d[7]*d[16]*d[25];
        coeffs[36] = d[4]*d[6]*d[13]*d[21] + d[5]*d[6]*d[16]*d[21] + d[4]*d[7]*d[13]*d[24] + d[5]*d[7]*d[16]*d[24];
        coeffs[37] = d[10]*d[20]*d[22] + d[11]*d[20]*d[25] + d[8]*d[14]*d[28] + d[9]*d[17]*d[28];
        coeffs[38] = -d[10]*d[18]*d[22] - d[11]*d[18]*d[25] - d[8]*d[12]*d[28] - d[9]*d[15]*d[28];
        coeffs[39] = -d[10]*d[19]*d[23] - d[11]*d[19]*d[26] - d[8]*d[13]*d[29] - d[9]*d[16]*d[29];
        coeffs[40] = d[10]*d[19]*d[22] + d[11]*d[19]*d[25] + d[8]*d[13]*d[28] + d[9]*d[16]*d[28];
        coeffs[41] = d[8]*d[10]*d[14]*d[22] + d[9]*d[10]*d[17]*d[22] + d[8]*d[11]*d[14]*d[25] + d[9]*d[11]*d[17]*d[25];  // NOLINT
        coeffs[42] = -d[8]*d[10]*d[12]*d[22] - d[9]*d[10]*d[15]*d[22] - d[8]*d[11]*d[12]*d[25]- d[9]*d[11]*d[15]*d[25];  // NOLINT
        coeffs[43] = -d[8]*d[10]*d[13]*d[23] - d[9]*d[10]*d[16]*d[23] - d[8]*d[11]*d[13]*d[26] - d[9]*d[11]*d[16]*d[26];  // NOLINT
        coeffs[44] = d[8]*d[10]*d[13]*d[22] + d[9]*d[10]*d[16]*d[22] + d[8]*d[11]*d[13]*d[25] + d[9]*d[11]*d[16]*d[25];

        // Setup elimination template
        static const int coeffs0_ind[] = { 2,3,2,3,2,3,14,3,14,3,6,7,3,2,23,24,37,2,2,3,7,17,14,3,24,31,38,3,3,14,18,15,32,15,8,19,16,4,25,33,39,4,4,16,9,5,26,40,5,5,10,11,7,6,27,28,41,23,37,24,11,20,17,7,28,34,42,24,38,31,21,18,35,32 };  // NOLINT
        static const int coeffs1_ind[] = { 11,10,27,41,28,20,11,28,42,34,21,35,22,12,29,43,36,12,22,19,8,29,36,43,25,39,33,13,30,44,13,9,30,44,26,40 };  // NOLINT
        static const int C0_ind[] = { 0,1,4,5,6,10,11,14,15,16,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,41,42,45,49,50,51,52,53,54,55,56,57,58,59,60,63,64,66,67,68,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,91,92,95,99 };  // NOLINT
        static const int C1_ind[] = { 2,3,7,8,9,12,13,17,18,19,22,29,32,33,37,38,39,40,41,42,43,44,45,46,47,48,49,53,57,58,60,63,64,66,67,68 };  // NOLINT

        Eigen::Matrix<double, 10, 10> C0; C0.setZero();
        Eigen::Matrix<double, 10, 7> C1; C1.setZero();
        for (int i = 0; i < 74; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
        for (int i = 0; i < 36; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); }

        Eigen::Matrix<double, 10, 7> C12 = C0.partialPivLu().solve(C1);



        // Setup action matrix
        Eigen::Matrix<double, 12, 7> RR;
        RR << -C12.bottomRows(5), Eigen::Matrix<double, 7, 7>::Identity(7, 7);

        static const int AM_ind[] = { 2,3,4,9,0,11,1 };  // NOLINT
        Eigen::Matrix<double, 7, 7> AM;
        for (int i = 0; i < 7; i++) {
            AM.row(i) = RR.row(AM_ind[i]);
        }

        Eigen::Matrix<std::complex<double>, 6, 7> sols;
        sols.setZero();

        // Solve eigenvalue problem
        Eigen::EigenSolver<Eigen::Matrix<double, 7, 7> > es(AM);
        Eigen::ArrayXcd D = es.eigenvalues();
        Eigen::ArrayXXcd V = es.eigenvectors();

        // Normalize eigenvectors
        Eigen::ArrayXcd normalization = (V.row(0).array().square() + V.row(1).array().square()).sqrt();

        for (int i = 0; i < 7; i++) {
            V.col(i) /= normalization(i);
        }

        sols.row(0) = V.row(0).array();
        sols.row(1) = V.row(1).array();
        sols.row(2) = V.row(2).array();
        sols.row(3) = V.row(3).array();
        sols.row(4) = V.row(5).array();
        sols.row(5) = D.transpose().array();

        return sols;
    }
}  // namespace HomLib::ValtonenOrnhagArxiv2020A
