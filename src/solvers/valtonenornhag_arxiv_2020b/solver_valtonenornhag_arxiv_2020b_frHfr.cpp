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

#include "solver_valtonenornhag_arxiv_2020b_frHfr.hpp"
#include <Eigen/Dense>

namespace HomLib::ValtonenOrnhagArxiv2020B {
    Eigen::MatrixXcd solver_frHfr(const Eigen::VectorXd& data) {
        // Compute coefficients
        const double* d = data.data();
        Eigen::VectorXd coeffs(64);
        double t1 = d[2]*d[36];
        double t2 = std::pow(d[20], 2);
        double t3 = std::pow(d[21], 2);
        double t4 = t2+t3;
        double t5 = std::pow(d[18], 2);
        double t6 = std::pow(d[19], 2);
        double t7 = t5+t6;
        double t8 = t7*t4;
        double t10 = t5+t6+t2+t3;
        double t12 = t4*d[5];
        double t16 = d[11]*d[36];
        double t25 = d[18]*d[34]+d[19]*d[35];
        double t37 = t12+d[2];
        double t44 = d[8]*d[36];
        double t50 = d[36]*d[17];
        double t51 = d[30]*d[31];
        double t52 = 2*t51;
        double t53 = d[32]*d[33];
        double t54 = 2*t53;
        double t55 = t50+t52-t54;
        double t56 = t5*t55;
        double t60 = d[11]*d[35];
        double t73 = d[21]*d[36];
        double t77 = t2*d[14]+t3*d[14]+d[11];
        double t84 = t53-0.5*t50-t51;
        double t95 = std::pow(d[30], 2);
        double t96 = std::pow(d[31], 2);
        double t97 = std::pow(d[32], 2);
        double t98 = std::pow(d[33], 2);
        double t99 = d[17]*d[35]-t95+t96-t97+t98;
        double t100 = d[19]*t99;
        double t106 = d[17]*d[34]-2.0*d[30]*d[33]-2*d[31]*d[32];
        double t107 = t106*d[18];
        double t108 = t100+t107;
        double t110 = t95-t96-t97+t98;
        double t129 = std::pow(d[24], 2);
        double t130 = std::pow(d[25], 2);
        double t131 = t129+t130;
        double t132 = std::pow(d[22], 2);
        double t133 = std::pow(d[23], 2);
        double t134 = t132+t133;
        double t135 = t134*t131;
        double t137 = t132+t133+t129+t130;
        double t139 = t131*d[5];
        double t151 = d[22]*d[34]+d[23]*d[35];
        double t161 = t139+d[2];
        double t173 = t132*t55;
        double t188 = d[25]*d[36];
        double t192 = t129*d[14]+t130*d[14]+d[11];
        double t207 = d[23]*t99;
        double t208 = t106*d[22];
        double t209 = t207+t208;
        coeffs[0] = t8*t1;
        coeffs[1] = (t10*d[2]+t12*t7)*d[36];
        coeffs[2] = t8*t16;
        coeffs[3] = (t10*d[5]+d[2])*d[36];
        coeffs[4] = t4*(t7*d[8]*d[36]+t25*d[2]);
        coeffs[5] = (t4*t7*d[14]+t10*d[11])*d[36];
        coeffs[6] = d[5]*d[36];
        coeffs[7] = t37*d[18]*d[34]+t37*d[19]*d[35]+t5*d[8]*d[36]+t6*d[8]*d[36]+t4*t44;
        coeffs[8] = (t10*d[14]+d[11])*d[36];
        coeffs[9] = t4*(t56+d[11]*d[18]*d[34]+d[19]*(t55*d[19]+t60));
        coeffs[10] = d[5]*d[18]*d[34]+d[5]*d[19]*d[35]+t44;
        coeffs[11] = d[14]*d[36];
        coeffs[12] = t25*t4*d[8];
        coeffs[13] = t7*t73;
        coeffs[14] = t77*d[18]*d[34]+t77*d[19]*d[35]-2.0*t4*t84+t55*t6+t56;
        coeffs[15] = t25*d[8];
        coeffs[16] = t73;
        coeffs[17] = d[14]*d[18]*d[34]+d[14]*d[19]*d[35]+t50+t52-t54;
        coeffs[18] = t110*t7*d[21]+t108*t2+t108*t3;
        coeffs[19] = t25*d[21];
        coeffs[20] = t110*d[21]+t100+t107;
        coeffs[21] = -2.0*((-d[18]*d[31]-d[19]*d[32])*d[33]+(d[18]*d[32]-d[19]*d[31])*d[30])*d[21];
        coeffs[22] = t135*t1;
        coeffs[23] = (t134*t139+t137*d[2])*d[36];
        coeffs[24] = t135*t16;
        coeffs[25] = (t137*d[5]+d[2])*d[36];
        coeffs[26] = t131*(t134*d[8]*d[36]+t151*d[2]);
        coeffs[27] = (t131*t134*d[14]+t137*d[11])*d[36];
        coeffs[28] = t133*d[8]*d[36]+t161*d[22]*d[34]+t161*d[23]*d[35]+t131*t44+t132*t44;
        coeffs[29] = (t137*d[14]+d[11])*d[36];
        coeffs[30] = t131*(t173+d[11]*d[34]*d[22]+d[23]*(t55*d[23]+t60));
        coeffs[31] = d[5]*d[22]*d[34]+d[5]*d[23]*d[35]+t44;
        coeffs[32] = t151*t131*d[8];
        coeffs[33] = t134*t188;
        coeffs[34] = t192*d[22]*d[34]+t192*d[23]*d[35]-2.0*t131*t84+t133*t55+t173;
        coeffs[35] = t151*d[8];
        coeffs[36] = t188;
        coeffs[37] = d[14]*d[22]*d[34]+d[14]*d[23]*d[35]+t50+t52-t54;
        coeffs[38] = t110*t134*d[25]+t129*t209+t130*t209;
        coeffs[39] = t151*d[25];
        coeffs[40] = t110*d[25]+t207+t208;
        coeffs[41] = -2.0*((-d[22]*d[31]-d[23]*d[32])*d[33]+(d[22]*d[32]-d[23]*d[31])*d[30])*d[25];
        coeffs[42] = d[1];
        coeffs[43] = -d[0]+d[4];
        coeffs[44] = d[10];
        coeffs[45] = -d[3];
        coeffs[46] = d[7];
        coeffs[47] = -d[9]+d[13];
        coeffs[48] = -d[6];
        coeffs[49] = -d[12];
        coeffs[50] = d[16];
        coeffs[51] = -d[15];
        coeffs[52] = d[2];
        coeffs[53] = d[5];
        coeffs[54] = -d[1];
        coeffs[55] = d[11];
        coeffs[56] = -d[4]+d[8];
        coeffs[57] = d[14];
        coeffs[58] = -d[10];
        coeffs[59] = -d[7];
        coeffs[60] = -d[13]+d[17];
        coeffs[61] = -d[16];
        coeffs[62] = -d[0]+d[8];
        coeffs[63] = -d[9]+d[17];

        // Setup elimination template
        static const int coeffs0_ind[] = { 0,22,52,42,52,1,22,23,0,52,53,43,42,52,53,3,23,25,1,52,53,4,26,0,42,52,62,22,45,43,52,53,6,25,6,3,52,53,46,42,52,62,7,26,28,1,4,43,42,53,52,62,45,23,45,53,48,46,54,43,62,52,53,45,10,28,31,3,7,45,52,54,43,62,53,45,25,12,32,4,46,42,52,62,48,26,13,50,33,44,55,63,48,56,45,45,53,31,6,10,53,56,45,45,6,15,32,35,7,12,48,54,46,43,53,45,62,48,52,28,16,51,33,50,36,58,13,47,63,55,57,49,18,38,9,50,44,55,63,51,30,35,10,15,56,59,48,45,48,45,53,31,36,51,60,16,49,49,57,12,46,62,48,54,32,19,39,13,50,58,63,33,51,20,38,40,14,18,51,58,50,47,57,49,63,51,55,34,40,17,20,60,61,51,49,51,49,57,37 };  // NOLINT
        static const int coeffs1_ind[] = { 21,41,18,50,63,51,58,38,41,20,21,61,51,49,51,60,40,21,51,61,41 };  // NOLINT
        static const int C0_ind[] = { 0,4,22,27,51,52,54,56,59,73,74,79,81,88,103,104,106,108,111,122,125,130,134,135,139,149,152,154,157,159,162,166,182,184,186,189,195,200,209,216,223,233,234,236,238,239,241,243,248,253,254,255,256,258,263,266,287,289,292,294,296,297,301,311,312,314,316,317,319,321,324,325,326,330,332,333,336,338,342,343,347,354,355,357,360,362,364,365,368,372,379,389,393,396,398,400,401,418,421,423,428,429,430,434,440,442,444,446,447,449,451,454,456,458,459,461,462,463,465,466,468,469,470,471,472,474,475,476,478,479,483,493,494,498,499,503,510,511,513,516,518,522,525,527,532,533,534,536,538,540,543,544,548,549,552,553,554,556,557,577,588,589,591,595,596,598,602,603,606,609,613,622,623,624,626,628,629,631,633,636,638,640,641,643,644,645,647,648,652,655,657,662,663,664,666,668,670,673,674 };  // NOLINT
        static const int C1_ind[] = { 0,4,5,16,17,19,23,24,28,31,33,38,42,43,46,49,50,57,69,75,76 };  // NOLINT

        Eigen::Matrix<double, 26, 26> C0; C0.setZero();
        Eigen::Matrix<double, 26, 3> C1; C1.setZero();
        for (int i = 0; i < 199; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
        for (int i = 0; i < 21; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); }

        Eigen::Matrix<double, 26, 3> C12 = C0.partialPivLu().solve(C1);

        // Setup action matrix
        Eigen::Matrix<double, 5, 3> RR;
        RR << -C12.bottomRows(2), Eigen::Matrix<double, 3, 3>::Identity(3, 3);

        static const int AM_ind[] = { 0,1,3 };  // NOLINT
        Eigen::Matrix<double, 3, 3> AM;
        for (int i = 0; i < 3; i++) {
            AM.row(i) = RR.row(AM_ind[i]);
        }

        // TODO(marcusvaltonen): Return only two rows.
        Eigen::Matrix<std::complex<double>, 4, 3> sols;
        sols.setZero();

        // Solve eigenvalue problem
        Eigen::EigenSolver<Eigen::Matrix<double, 3, 3> > es(AM);
        Eigen::ArrayXcd D = es.eigenvalues();
        Eigen::ArrayXXcd V = es.eigenvectors();

        V = (V / V.row(2).array().replicate(3, 1)).eval();

        sols.row(2) = D.transpose().array();
        sols.row(3) = V.row(0).array() / (sols.row(2).array());

        return sols;
    }
}  // namespace HomLib::ValtonenOrnhagArxiv2020B
