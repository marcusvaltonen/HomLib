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
#include "PoseLib/misc/sturm.h"

namespace HomLib {
namespace KukelovaCVPR2015 {
    int solver_kukelova_cvpr_2015(const Eigen::VectorXd& data, Eigen::Matrix<double, 3, 5>& sols) {
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
        RR.topRows(4) = -C12.bottomRows(4);
        RR.bottomRows(5) = Eigen::Matrix<double, 5, 5>::Identity(5, 5);

        static const int AM_ind[] = { 5,0,1,2,3 };  // NOLINT
        Eigen::Matrix<double, 5, 5> AM;
        for (int i = 0; i < 5; i++) {
            AM.row(i) = RR.row(AM_ind[i]);
        }

        // Find real eigenvalues via characteristic polynomial + Sturm bisection. tol = 1e-14: with
        // the default 1e-10, Newton's method stops when |p(x)| < 1e-10, and the ill-conditioned roots
        // stay inaccurate.
        double p[6];
        Eigen::Matrix<double, 5, 5> AMp = AM;
        poselib::sturm::charpoly_danilevsky_piv(AMp, p);
        double roots[5];
        int nroots = poselib::sturm::bisect_sturm<5>(p, roots, 1e-14);

        // Close eigenvalues are ill-conditioned roots of the characteristic polynomial, which only
        // determines them roughly: then Eigen's eigensolver, as before (about 1 problem in 100).
        // The roots are in increasing order.
        for (int i = 0; i + 1 < nroots; i++) {
            if (roots[i + 1] - roots[i] < 1e-3 * (1.0 + std::abs(roots[i]))) {
                Eigen::EigenSolver<Eigen::Matrix<double, 5, 5>> es(AM);
                if (es.info() != Eigen::Success) {
                    return 0;
                }
                int n = 0;
                for (int j = 0; j < 5; j++) {
                    const std::complex<double> v0 = es.eigenvectors()(0, j);
                    const Eigen::Vector3cd sol(es.eigenvalues()(j), es.eigenvectors()(2, j) / v0,
                                               es.eigenvectors()(3, j) / v0);
                    // The complex solutions are removed, with the same threshold as before
                    if (sol.imag().squaredNorm() <= 1e-5) {
                        sols.col(n++) = sol.real();
                    }
                }
                return n;
            }
        }

        // Fast eigenvector extraction for each real eigenvalue.
        // AM.row(0) = [0,1,0,0,0], so v[1] = λ·v[0]; with v[0] = 1: v[1] = λ.
        // Rows 2–4 of (AM − λI)v = 0 give a 3×3 system for v[2], v[3], v[4].
        for (int i = 0; i < nroots; i++) {
            const double lam = roots[i];
            sols(0, i) = lam;

            Eigen::Matrix3d B;
            Eigen::Vector3d rhs;
            for (int r = 0; r < 3; r++) {
                B(r, 0) = AM(r + 2, 2);
                B(r, 1) = AM(r + 2, 3);
                B(r, 2) = AM(r + 2, 4);
                B(r, r) -= lam;
                rhs(r) = -(AM(r + 2, 0) + AM(r + 2, 1) * lam);
            }
            Eigen::Vector3d v = B.partialPivLu().solve(rhs);

            // lam is accurate only up to the conditioning of the characteristic polynomial: up to 2
            // steps of Newton's method on (AM - lam I) w = 0, w = (1, lam, v), rows 1-4 (row 0 holds
            // by construction), with unknowns lam and v, each kept only if it reduces the residual
            double l = lam;
            Eigen::Vector4d F;
            Eigen::Matrix4d J;
            const auto residual = [&AM](double ll, const Eigen::Vector3d &vv, Eigen::Vector4d &FF,
                                        Eigen::Matrix4d &JJ) {
                const Eigen::Matrix<double, 5, 1> w(1.0, ll, vv(0), vv(1), vv(2));
                for (int r = 1; r < 5; r++) {
                    FF(r - 1) = AM.row(r).dot(w) - ll * w(r);
                    JJ(r - 1, 0) = AM(r, 1) - w(r) - (r == 1 ? ll : 0.0);
                    for (int k = 2; k < 5; k++)
                        JJ(r - 1, k - 1) = AM(r, k) - (r == k ? ll : 0.0);
                }
            };
            residual(l, v, F, J);
            for (int it = 0; it < 2; it++) {
                const Eigen::Vector4d du = J.partialPivLu().solve(-F);
                const double l1 = l + du(0);
                const Eigen::Vector3d v1 = v + du.tail<3>();
                Eigen::Vector4d F1;
                Eigen::Matrix4d J1;
                residual(l1, v1, F1, J1);
                if (!(F1.norm() < F.norm())) {
                    break;
                }
                l = l1;
                v = v1;
                F = F1;
                J = J1;
            }
            sols(0, i) = l;
            sols(1, i) = v(0);
            sols(2, i) = v(1);
        }

        return nroots;
    }
}  // namespace KukelovaCVPR2015
}  // namespace HomLib
