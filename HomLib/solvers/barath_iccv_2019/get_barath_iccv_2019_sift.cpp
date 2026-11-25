// Copyright (C) 2019 Czech Technical University.
// All rights reserved.
//
// Port of HomographyTwoSIFTSolver from pygcransac (Barath & Kukelova, ICCV 2019)
// to HomLib's Eigen-based API.  The algorithm is kept exactly as in the original.

#include <vector>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include "posedata.hpp"
#include "get_barath_iccv_2019.hpp"

namespace HomLib {
namespace BarathICCV2019 {

std::vector<HomLib::PoseData> get_sift(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    const std::vector<Eigen::Vector2d> &ori,
    const std::vector<Eigen::Vector2d> &scale
) {
    std::vector<HomLib::PoseData> output;
    if (x.size() < 2)
        return output;

    constexpr int kRowNumber = 6;
    Eigen::Matrix<double, 6, 9> coefficients;

    double tmp_data[20];
    int data_idx = 0;
    int rowIdx = 0;

    for (int i = 0; i < 2; ++i) {
        const double s1 = std::sin(ori[i](0)), c1 = std::cos(ori[i](0));
        const double s2 = std::sin(ori[i](1)), c2 = std::cos(ori[i](1));

        const double u1 = x[i](0), v1 = x[i](1);
        const double u2 = y[i](0), v2 = y[i](1);
        const double q_ratio = scale[i](1) / scale[i](0);
        const double q1 = 1.0;
        const double q2 = q_ratio * q_ratio;

        coefficients.row(rowIdx++) << 0, 0, 0, u1, v1, 1, -u1 * v2, -v1 * v2, -v2;
        coefficients.row(rowIdx++) << u1, v1, 1, 0, 0, 0, -u1 * u2, -v1 * u2, -u2;

        tmp_data[data_idx++] = u1;
        tmp_data[data_idx++] = v1;
        tmp_data[data_idx++] = u2;
        tmp_data[data_idx++] = v2;
        tmp_data[data_idx++] = s1;
        tmp_data[data_idx++] = c1;
        tmp_data[data_idx++] = s2;
        tmp_data[data_idx++] = c2;
        tmp_data[data_idx++] = q1;
        tmp_data[data_idx++] = q2;

        const double c1s2 = c1 * s2;
        const double c1c2 = c1 * c2;
        const double s1s2 = s1 * s2;
        const double s1c2 = s1 * c2;

        coefficients.row(rowIdx++) << -c1s2, -s1s2, 0, c1c2, s1c2, 0,
            u2 * c1s2 - v2 * c1c2, u2 * s1s2 - v2 * s1c2, 0;
    }

    Eigen::JacobiSVD<Eigen::Matrix<double, 6, 9>> svd(
        coefficients, Eigen::ComputeFullV);
    Eigen::Matrix<double, 9, 3> null_space = svd.matrixV().rightCols<3>();

    const double n11 = null_space(0, 0), n12 = null_space(1, 0), n13 = null_space(2, 0);
    const double n14 = null_space(3, 0), n15 = null_space(4, 0), n16 = null_space(5, 0);
    const double n17 = null_space(6, 0), n18 = null_space(7, 0), n19 = null_space(8, 0);
    const double n21 = null_space(0, 1), n22 = null_space(1, 1), n23 = null_space(2, 1);
    const double n24 = null_space(3, 1), n25 = null_space(4, 1), n26 = null_space(5, 1);
    const double n27 = null_space(6, 1), n28 = null_space(7, 1), n29 = null_space(8, 1);
    const double n31 = null_space(0, 2), n32 = null_space(1, 2), n33 = null_space(2, 2);
    const double n34 = null_space(3, 2), n35 = null_space(4, 2), n36 = null_space(5, 2);
    const double n37 = null_space(6, 2), n38 = null_space(7, 2), n39 = null_space(8, 2);

    const double d1 = tmp_data[0], d2 = tmp_data[1], d3 = tmp_data[2], d4 = tmp_data[3];
    const double d9 = tmp_data[8], d10 = tmp_data[9];
    const double d11 = tmp_data[10], d12 = tmp_data[11], d13 = tmp_data[12], d14 = tmp_data[13];
    const double d19 = tmp_data[18], d20 = tmp_data[19];

    const double d1d10 = d1 * d10, d2d10 = d2 * d10;
    const double d1d1d10 = d1 * d1d10, d1d2d10 = d1 * d2d10, d2d2d10 = d2 * d2d10;
    const double d3d9 = d3 * d9, d4d9 = d4 * d9;
    const double d11d20 = d11 * d20, d12d20 = d12 * d20;
    const double d11d11d20 = d11 * d11d20, d11d12d20 = d11 * d12d20, d12d12d20 = d12 * d12d20;
    const double d13d19 = d13 * d19, d14d19 = d14 * d19;

    const double n15n17 = n15 * n17, n14n18 = n14 * n18;
    const double n12n17 = n12 * n17, n11n18 = n11 * n18;
    const double n17n17 = n17 * n17, n17n18 = n17 * n18, n18n18 = n18 * n18;
    const double n18n24 = n18 * n24, n17n25 = n17 * n25;
    const double n15n27 = n15 * n27, n14n28 = n14 * n28;
    const double n18n21 = n18 * n21, n17n22 = n17 * n22;
    const double n25n27 = n25 * n27, n24n28 = n24 * n28;
    const double n22n27 = n22 * n27, n21n28 = n21 * n28;
    const double n27n27 = n27 * n27, n27n28 = n27 * n28, n28n28 = n28 * n28;
    const double n18n34 = n18 * n34, n17n35 = n17 * n35;
    const double n15n37 = n15 * n37, n14n38 = n14 * n38;
    const double n18n31 = n18 * n31, n17n32 = n17 * n32;
    const double n28n34 = n28 * n34, n27n35 = n27 * n35;
    const double n25n37 = n25 * n37, n24n38 = n24 * n38;
    const double n28n31 = n28 * n31, n27n32 = n27 * n32;
    const double n35n37 = n35 * n37, n34n38 = n34 * n38;
    const double n32n37 = n32 * n37, n31n38 = n31 * n38;
    const double n37n37 = n37 * n37, n37n38 = n37 * n38, n38n38 = n38 * n38;
    const double n12n14 = n12 * n14, n11n15 = n11 * n15;
    const double n17n19 = n17 * n19, n18n19 = n18 * n19, n19n19 = n19 * n19;
    const double n12n27 = n12 * n27, n11n28 = n11 * n28;
    const double n17n27 = n17 * n27, n18n27 = n18 * n27;
    const double n17n28 = n17 * n28, n18n28 = n18 * n28;
    const double n15n21 = n15 * n21, n14n22 = n14 * n22;
    const double n12n24 = n12 * n24, n11n25 = n11 * n25;
    const double n19n27 = n19 * n27, n17n29 = n17 * n29;
    const double n19n28 = n19 * n28, n18n29 = n18 * n29;
    const double n19n29 = n19 * n29;
    const double n22n24 = n22 * n24, n21n25 = n21 * n25;
    const double n27n29 = n27 * n29, n28n29 = n28 * n29, n29n29 = n29 * n29;
    const double n12n37 = n12 * n37, n11n38 = n11 * n38;
    const double n17n37 = n17 * n37, n18n37 = n18 * n37;
    const double n17n38 = n17 * n38, n18n38 = n18 * n38;
    const double n15n31 = n15 * n31, n14n32 = n14 * n32;
    const double n12n34 = n12 * n34, n11n35 = n11 * n35;
    const double n19n37 = n19 * n37, n17n39 = n17 * n39;
    const double n19n38 = n19 * n38, n18n39 = n18 * n39;
    const double n19n39 = n19 * n39;
    const double n22n37 = n22 * n37, n21n38 = n21 * n38;
    const double n27n37 = n27 * n37, n28n37 = n28 * n37;
    const double n27n38 = n27 * n38, n28n38 = n28 * n38;
    const double n25n31 = n25 * n31, n24n32 = n24 * n32;
    const double n22n34 = n22 * n34, n21n35 = n21 * n35;
    const double n29n37 = n29 * n37, n27n39 = n27 * n39;
    const double n29n38 = n29 * n38, n28n39 = n28 * n39;
    const double n29n39 = n29 * n39;
    const double n32n34 = n32 * n34, n31n35 = n31 * n35;
    const double n37n39 = n37 * n39, n38n39 = n38 * n39, n39n39 = n39 * n39;

    // Constraint coefficients (from correspondence 1)
    const double coeffs0 = n15n17 * d3d9 - n14n18 * d3d9 - n12n17 * d4d9 + n11n18 * d4d9 + n17n17 * d1d1d10 + 2 * n17n18 * d1d2d10 + n18n18 * d2d2d10 + n12n14 * d9 - n11n15 * d9 + 2 * n17n19 * d1d10 + 2 * n18n19 * d2d10 + n19n19 * d10;
    const double coeffs1 = -n18n24 * d3d9 + n17n25 * d3d9 + n15n27 * d3d9 - n14n28 * d3d9 + n18n21 * d4d9 - n17n22 * d4d9 - n12n27 * d4d9 + n11n28 * d4d9 + 2 * n17n27 * d1d1d10 + 2 * n18n27 * d1d2d10 + 2 * n17n28 * d1d2d10 + 2 * n18n28 * d2d2d10 - n15n21 * d9 + n14n22 * d9 + n12n24 * d9 - n11n25 * d9 + 2 * n19n27 * d1d10 + 2 * n17n29 * d1d10 + 2 * n19n28 * d2d10 + 2 * n18n29 * d2d10 + 2 * n19n29 * d10;
    const double coeffs2 = n25n27 * d3d9 - n24n28 * d3d9 - n22n27 * d4d9 + n21n28 * d4d9 + n27n27 * d1d1d10 + 2 * n27n28 * d1d2d10 + n28n28 * d2d2d10 + n22n24 * d9 - n21n25 * d9 + 2 * n27n29 * d1d10 + 2 * n28n29 * d2d10 + n29n29 * d10;
    const double coeffs3 = -n18n34 * d3d9 + n17n35 * d3d9 + n15n37 * d3d9 - n14n38 * d3d9 + n18n31 * d4d9 - n17n32 * d4d9 - n12n37 * d4d9 + n11n38 * d4d9 + 2 * n17n37 * d1d1d10 + 2 * n18n37 * d1d2d10 + 2 * n17n38 * d1d2d10 + 2 * n18n38 * d2d2d10 - n15n31 * d9 + n14n32 * d9 + n12n34 * d9 - n11n35 * d9 + 2 * n19n37 * d1d10 + 2 * n17n39 * d1d10 + 2 * n19n38 * d2d10 + 2 * n18n39 * d2d10 + 2 * n19n39 * d10;
    const double coeffs4 = -n28n34 * d3d9 + n27n35 * d3d9 + n25n37 * d3d9 - n24n38 * d3d9 + n28n31 * d4d9 - n27n32 * d4d9 - n22n37 * d4d9 + n21n38 * d4d9 + 2 * n27n37 * d1d1d10 + 2 * n28n37 * d1d2d10 + 2 * n27n38 * d1d2d10 + 2 * n28n38 * d2d2d10 - n25n31 * d9 + n24n32 * d9 + n22n34 * d9 - n21n35 * d9 + 2 * n29n37 * d1d10 + 2 * n27n39 * d1d10 + 2 * n29n38 * d2d10 + 2 * n28n39 * d2d10 + 2 * n29n39 * d10;
    const double coeffs5 = n35n37 * d3d9 - n34n38 * d3d9 - n32n37 * d4d9 + n31n38 * d4d9 + n37n37 * d1d1d10 + 2 * n37n38 * d1d2d10 + n38n38 * d2d2d10 + n32n34 * d9 - n31n35 * d9 + 2 * n37n39 * d1d10 + 2 * n38n39 * d2d10 + n39n39 * d10;

    // Constraint coefficients (from correspondence 2)
    const double coeffs6 = n15n17 * d13d19 - n14n18 * d13d19 - n12n17 * d14d19 + n11n18 * d14d19 + n17n17 * d11d11d20 + 2 * n17n18 * d11d12d20 + n18n18 * d12d12d20 + n12n14 * d19 - n11n15 * d19 + 2 * n17n19 * d11d20 + 2 * n18n19 * d12d20 + n19n19 * d20;
    const double coeffs7 = -n18n24 * d13d19 + n17n25 * d13d19 + n15n27 * d13d19 - n14n28 * d13d19 + n18n21 * d14d19 - n17n22 * d14d19 - n12n27 * d14d19 + n11n28 * d14d19 + 2 * n17n27 * d11d11d20 + 2 * n18n27 * d11d12d20 + 2 * n17n28 * d11d12d20 + 2 * n18n28 * d12d12d20 - n15n21 * d19 + n14n22 * d19 + n12n24 * d19 - n11n25 * d19 + 2 * n19n27 * d11d20 + 2 * n17n29 * d11d20 + 2 * n19n28 * d12d20 + 2 * n18n29 * d12d20 + 2 * n19n29 * d20;
    const double coeffs8 = n25n27 * d13d19 - n24n28 * d13d19 - n22n27 * d14d19 + n21n28 * d14d19 + n27n27 * d11d11d20 + 2 * n27n28 * d11d12d20 + n28n28 * d12d12d20 + n22n24 * d19 - n21n25 * d19 + 2 * n27n29 * d11d20 + 2 * n28n29 * d12d20 + n29n29 * d20;
    const double coeffs9 = -n18n34 * d13d19 + n17n35 * d13d19 + n15n37 * d13d19 - n14n38 * d13d19 + n18n31 * d14d19 - n17n32 * d14d19 - n12n37 * d14d19 + n11n38 * d14d19 + 2 * n17n37 * d11d11d20 + 2 * n18n37 * d11d12d20 + 2 * n17n38 * d11d12d20 + 2 * n18n38 * d12d12d20 - n15n31 * d19 + n14n32 * d19 + n12n34 * d19 - n11n35 * d19 + 2 * n19n37 * d11d20 + 2 * n17n39 * d11d20 + 2 * n19n38 * d12d20 + 2 * n18n39 * d12d20 + 2 * n19n39 * d20;
    const double coeffs10 = -n28n34 * d13d19 + n27n35 * d13d19 + n25n37 * d13d19 - n24n38 * d13d19 + n28n31 * d14d19 - n27n32 * d14d19 - n22n37 * d14d19 + n21n38 * d14d19 + 2 * n27n37 * d11d11d20 + 2 * n28n37 * d11d12d20 + 2 * n27n38 * d11d12d20 + 2 * n28n38 * d12d12d20 - n25n31 * d19 + n24n32 * d19 + n22n34 * d19 - n21n35 * d19 + 2 * n29n37 * d11d20 + 2 * n27n39 * d11d20 + 2 * n29n38 * d12d20 + 2 * n28n39 * d12d20 + 2 * n29n39 * d20;
    const double coeffs11 = n35n37 * d13d19 - n34n38 * d13d19 - n32n37 * d14d19 + n31n38 * d14d19 + n37n37 * d11d11d20 + 2 * n37n38 * d11d12d20 + n38n38 * d12d12d20 + n32n34 * d19 - n31n35 * d19 + 2 * n37n39 * d11d20 + 2 * n38n39 * d12d20 + n39n39 * d20;

    // Fill the 10x6 constraint matrix C (row-major flat array)
    double C[60];
    std::fill(C, C + 60, 0.0);

    C[0] = coeffs0;  C[1] = coeffs1;  C[2] = coeffs2;  C[3] = coeffs3;  C[4] = coeffs4;
    C[7] = coeffs5;
    C[11] = coeffs0;
    C[12] = coeffs1;  C[14] = coeffs3;  C[15] = coeffs2;
    C[18] = coeffs5;  C[19] = coeffs4;
    C[21] = coeffs6;  C[22] = coeffs7;  C[24] = coeffs9;  C[25] = coeffs8;
    C[28] = coeffs11; C[29] = coeffs10;
    C[33] = coeffs6;  C[34] = coeffs7;
    C[36] = coeffs11; C[37] = coeffs9;  C[38] = coeffs10; C[39] = coeffs8;
    C[43] = coeffs0;  C[44] = coeffs1;  C[46] = coeffs5;  C[47] = coeffs3;
    C[48] = coeffs4;  C[49] = coeffs2;  C[50] = coeffs6;  C[51] = coeffs7;  C[52] = coeffs8;  C[53] = coeffs9;
    C[54] = coeffs10; C[57] = coeffs11;

    // Extract C0 (6x6) and C1 (6x4) from C using the same byte-offset pattern as the original
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> C0;
    Eigen::Matrix<double, 6, 4, Eigen::RowMajor> C1;

    std::memcpy(C0.data(),      C,      6 * sizeof(double));
    std::memcpy(C0.data() + 6,  C + 10, 6 * sizeof(double));
    std::memcpy(C0.data() + 12, C + 20, 6 * sizeof(double));
    std::memcpy(C0.data() + 18, C + 30, 6 * sizeof(double));
    std::memcpy(C0.data() + 24, C + 40, 6 * sizeof(double));
    std::memcpy(C0.data() + 30, C + 50, 6 * sizeof(double));

    std::memcpy(C1.data(),      C + 6,  4 * sizeof(double));
    std::memcpy(C1.data() + 4,  C + 16, 4 * sizeof(double));
    std::memcpy(C1.data() + 8,  C + 26, 4 * sizeof(double));
    std::memcpy(C1.data() + 12, C + 36, 4 * sizeof(double));
    std::memcpy(C1.data() + 16, C + 46, 4 * sizeof(double));
    std::memcpy(C1.data() + 20, C + 56, 4 * sizeof(double));

    // Solve C0 * C12 = C1
    Eigen::Matrix<double, 6, 4> C12 = C0.fullPivLu().solve(C1);

    // Build action matrix
    // RR = [-C12 bottom 2 rows; I_4]
    Eigen::Matrix<double, 6, 4, Eigen::RowMajor> RR;
    RR.setZero();
    RR.row(0) = -C12.row(4);
    RR.row(1) = -C12.row(5);
    RR(2, 0) = 1.0;
    RR(3, 1) = 1.0;
    RR(4, 2) = 1.0;
    RR(5, 3) = 1.0;

    // AM = RR(AM_ind, :) where AM_ind = {4, 0, 5, 1}
    Eigen::Matrix4d AM;
    AM.row(0) = RR.row(4);
    AM.row(1) = RR.row(0);
    AM.row(2) = RR.row(5);
    AM.row(3) = RR.row(1);

    // Eigenvalue decomposition
    Eigen::EigenSolver<Eigen::Matrix4d> es(AM);
    if (es.info() != Eigen::Success)
        return output;

    Eigen::Vector4cd D = es.eigenvalues();
    Eigen::Matrix4cd V = es.eigenvectors();

    // Normalize eigenvectors: V = V ./ V(0,:)
    for (int i = 0; i < 4; ++i) {
        if (std::abs(V(0, i)) < 1e-14)
            continue;
        V.col(i) /= V(0, i);
    }

    for (int i = 0; i < 4; ++i) {
        if (std::abs(D(i).imag()) > 1e-8 || std::abs(V(1, i).imag()) > 1e-8)
            continue;

        const double alpha = V(1, i).real();
        const double beta = D(i).real();

        Eigen::Matrix3d H;
        H(0, 0) = alpha * n11 + beta * n21 + n31;
        H(0, 1) = alpha * n12 + beta * n22 + n32;
        H(0, 2) = alpha * n13 + beta * n23 + n33;
        H(1, 0) = alpha * n14 + beta * n24 + n34;
        H(1, 1) = alpha * n15 + beta * n25 + n35;
        H(1, 2) = alpha * n16 + beta * n26 + n36;
        H(2, 0) = alpha * n17 + beta * n27 + n37;
        H(2, 1) = alpha * n18 + beta * n28 + n38;
        H(2, 2) = alpha * n19 + beta * n29 + n39;

        if (H(2, 2) != 0.0)
            H /= H(2, 2);
        if (!H.allFinite())
            continue;

        HomLib::PoseData pd;
        pd.homography = H;
        pd.distortion_parameter = 0.0;
        pd.distortion_parameter2 = 0.0;
        output.push_back(pd);
    }

    return output;
}

}  // namespace BarathICCV2019
}  // namespace HomLib
