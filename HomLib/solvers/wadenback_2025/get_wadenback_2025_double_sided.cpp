#include <Eigen/Dense>
#include <vector>

#include "PoseLib/misc/sturm.h"

#include "posedata.hpp"
#include "wadenback_common.hpp"
#include <iostream>

namespace HomLib {
namespace Wadenback2025 {
static inline void fast_eigenvector_solver(
    double *eigv,
    int neig,
    Eigen::Matrix<double, 9, 9> &AM,
    Eigen::Matrix<double, 2, 9> &sols
);
std::vector<HomLib::PoseData> get_double_sided(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y,
    bool extra_check
)
{
    Eigen::Matrix<double, 24, 1> data;
    data << HomLib::Wadenback2025::compute_coeffs(y), HomLib::Wadenback2025::compute_coeffs(x);  // NOTE: y first

	// Compute coefficients
    const double* d = data.data();
    Eigen::VectorXd coeffs(48);
    coeffs[0] = -d[2]*d[13] + d[1]*d[14];
    coeffs[1] = -d[5]*d[13] + d[4]*d[14];
    coeffs[2] = -d[2]*d[16] + d[1]*d[17];
    coeffs[3] = -d[8]*d[13] + d[7]*d[14];
    coeffs[4] = -d[5]*d[16] + d[4]*d[17];
    coeffs[5] = -d[2]*d[19] + d[1]*d[20];
    coeffs[6] = -d[11]*d[13] + d[10]*d[14];
    coeffs[7] = -d[8]*d[16] + d[7]*d[17];
    coeffs[8] = -d[5]*d[19] + d[4]*d[20];
    coeffs[9] = -d[2]*d[22] + d[1]*d[23];
    coeffs[10] = -d[11]*d[16] + d[10]*d[17];
    coeffs[11] = -d[8]*d[19] + d[7]*d[20];
    coeffs[12] = -d[5]*d[22] + d[4]*d[23];
    coeffs[13] = -d[11]*d[19] + d[10]*d[20];
    coeffs[14] = -d[8]*d[22] + d[7]*d[23];
    coeffs[15] = -d[11]*d[22] + d[10]*d[23];
    coeffs[16] = d[2]*d[12] - d[0]*d[14];
    coeffs[17] = d[5]*d[12] - d[3]*d[14];
    coeffs[18] = d[2]*d[15] - d[0]*d[17];
    coeffs[19] = d[8]*d[12] - d[6]*d[14];
    coeffs[20] = d[5]*d[15] - d[3]*d[17];
    coeffs[21] = d[2]*d[18] - d[0]*d[20];
    coeffs[22] = d[11]*d[12] - d[9]*d[14];
    coeffs[23] = d[8]*d[15] - d[6]*d[17];
    coeffs[24] = d[5]*d[18] - d[3]*d[20];
    coeffs[25] = d[2]*d[21] - d[0]*d[23];
    coeffs[26] = d[11]*d[15] - d[9]*d[17];
    coeffs[27] = d[8]*d[18] - d[6]*d[20];
    coeffs[28] = d[5]*d[21] - d[3]*d[23];
    coeffs[29] = d[11]*d[18] - d[9]*d[20];
    coeffs[30] = d[8]*d[21] - d[6]*d[23];
    coeffs[31] = d[11]*d[21] - d[9]*d[23];
    coeffs[32] = -d[1]*d[12] + d[0]*d[13];
    coeffs[33] = -d[4]*d[12] + d[3]*d[13];
    coeffs[34] = -d[1]*d[15] + d[0]*d[16];
    coeffs[35] = -d[7]*d[12] + d[6]*d[13];
    coeffs[36] = -d[4]*d[15] + d[3]*d[16];
    coeffs[37] = -d[1]*d[18] + d[0]*d[19];
    coeffs[38] = -d[10]*d[12] + d[9]*d[13];
    coeffs[39] = -d[7]*d[15] + d[6]*d[16];
    coeffs[40] = -d[4]*d[18] + d[3]*d[19];
    coeffs[41] = -d[1]*d[21] + d[0]*d[22];
    coeffs[42] = -d[10]*d[15] + d[9]*d[16];
    coeffs[43] = -d[7]*d[18] + d[6]*d[19];
    coeffs[44] = -d[4]*d[21] + d[3]*d[22];
    coeffs[45] = -d[10]*d[18] + d[9]*d[19];
    coeffs[46] = -d[7]*d[21] + d[6]*d[22];
    coeffs[47] = -d[10]*d[21] + d[9]*d[22];


    // Setup elimination template
    static const int coeffs0_ind[] = { 0,16,32,1,0,17,16,32,33,2,18,34,3,1,19,16,17,0,32,33,35,4,2,20,18,34,36,7,4,23,18,20,2,34,36,39,6,3,22,17,19,1,33,35,38,6,19,22,3,35,38,22,6,38 };
    static const int coeffs1_ind[] = { 10,7,26,20,23,4,36,39,42,10,23,26,7,39,42,13,11,29,24,27,8,40,43,45,26,10,42,13,27,29,11,43,45,15,14,31,28,30,12,44,46,47,29,13,45,15,30,31,14,46,47,31,15,47 };
    static const int C0_ind[] = { 0,2,8,9,10,11,13,16,17,18,20,26,27,28,29,30,31,32,33,34,35,36,37,38,40,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,64,66,67,68,69,70,75,77,78 } ;
    static const int C1_ind[] = { 0,1,2,3,4,5,6,7,8,10,12,13,14,15,16,18,19,20,21,22,23,24,25,26,30,32,33,37,39,40,41,42,43,45,46,47,48,49,50,51,52,53,57,59,60,64,66,67,68,69,70,75,77,78 };

    Eigen::Matrix<double,9,9> C0; C0.setZero();
    Eigen::Matrix<double,9,9> C1; C1.setZero();
    for (int i = 0; i < 54; i++) { C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); }
    for (int i = 0; i < 54; i++) { C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); } 

    //std::cout << "C0 =\n" << C0 << std::endl;
    //std::cout << "C1 =\n" << C1 << std::endl;

    //Eigen::Matrix<double,9,9> C12 = C0.partialPivLu().solve(C1);  // partialPivLu seems to be the fastest from the standard solvers
    //std::cout << "C12 =\n" << C12 << std::endl;
    // Eigen::Matrix<double,9,9> C12 = C0.householderQr().solve(C1); // TODO: We only need the three bottom rows... why compute all?
    
    // Schur complement trick to avoid computing unused lines
    Eigen::Matrix<double, 6, 6> A = C0.block(0,0,6,6);
    Eigen::Matrix<double, 6, 3> B = C0.block(0,6,6,3);
    Eigen::Matrix<double, 3, 6> C = C0.block(6,0,3,6);
    Eigen::Matrix<double, 3, 3> D = C0.block(6,6,3,3);
    
    Eigen::Matrix<double, 3, 3> M = D - C * A.partialPivLu().solve(B);
    Eigen::Matrix<double, 3, 9> N = C1.bottomRows(3) - C * A.partialPivLu().solve(C1.topRows(6));
    Eigen::Matrix<double, 3, 9> C12b = M.partialPivLu().solve(N);
    
    //std::cout << "C12b =\n" << C12b << std::endl;
    // Setup action matrix
    Eigen::Matrix<double, 12, 9> RR;
    //RR.topRows(3) = -C12.bottomRows(3);
    RR.topRows(3) = -C12b;
    RR.bottomRows(9) = Eigen::Matrix<double,9,9>::Identity(9, 9);

    static const int AM_ind[] = { 0,1,3,2,4,5,6,7,9 };
    Eigen::Matrix<double, 9, 9> AM;
    for (int i = 0; i < 9; i++) {
	    AM.row(i) = RR.row(AM_ind[i]);
    }

    Eigen::Matrix<double, 2, 9> sols;
    sols.setZero();

    // Solve eigenvalue problem
    double p[10];
    Eigen::Matrix<double, 9, 9> AMp = AM;

    // Find real eigenvalues (roots)
    poselib::sturm::charpoly_danilevsky_piv(AMp, p);
    double roots[9];
    int nroots = poselib::sturm::bisect_sturm<9>(p, roots);
    fast_eigenvector_solver(roots, nroots, AM, sols);


    // TODO: This static threshold does not work - find a way to dynamically adapt
    double THRESH = 1e-6;

    std::vector<HomLib::PoseData> output;
    
    for (int i = 0; i < nroots; i++) {
        double k1 = sols(0, i);
        double k2 = sols(1, i);

        if (k1 > 0.1 || k1 < -1)
            continue;
        if (k2 > 0.1 || k2 < -1)
            continue;

        if (extra_check && nroots > 1) {
            double res1a = compute_residual1(x, k1);
            if (std::abs(res1a) < THRESH) {
                continue;
            }
            double res1b = compute_residual1(y, k2);
            if (std::abs(res1b) < THRESH) {
                continue;
            }
            double res2a = compute_residual2(x, k1);
            if (std::abs(res2a) < THRESH) {
                continue;
            }

            double res2b = compute_residual2(y, k2);
            if (std::abs(res2b) < THRESH) {
                continue;
            }
	    }
        // Compute homography
        Eigen::Matrix3d H;

        std::vector<Eigen::Vector2d> x_undist;
        std::vector<Eigen::Vector2d> y_undist;
        for (int i = 0; i < 4; i++) {
	        Eigen::Vector2d tmp = x[i] / (1 + k1 * x[i].squaredNorm());
	        x_undist.push_back(tmp);
	        tmp = y[i] / (1 + k2 * y[i].squaredNorm());
	        y_undist.push_back(tmp);
        }
        H = HomLib::Wadenback2025::homography_4pt_guo(x_undist, y_undist);

        // Package
        HomLib::PoseData container;
        container.homography = H;
        container.distortion_parameter = k1;
        container.distortion_parameter2 = k2;
        output.push_back(container);
    }


	return output;
}
// Action =  x
// Quotient ring basis (V) = x^2*y^2,x^2*y,x*y^2,x^2,x*y,y^2,x,y,1,
// Available monomials (RR*V) = x^3*y^2,x^3*y,x^3,x^2*y^2,x^2*y,x*y^2,x^2,x*y,y^2,x,y,1,


static void fast_eigenvector_solver(
    double *eigv,
    int neig,
    Eigen::Matrix<double, 9, 9> &AM,
    Eigen::Matrix<double, 2, 9> &sols
) {
    static const int ind[] = { 0,1 };	
    // Truncated action matrix containing non-trivial rows
    Eigen::Matrix<double, 2, 9> AMs;
    double zi[3];

    for (int i = 0; i < 2; i++)	{
	    AMs.row(i) = AM.row(ind[i]);
    }
    for (int i = 0; i < neig; i++) {
        zi[0] = eigv[i];
        for (int j = 1; j < 3; j++)
        {
            zi[j] = zi[j - 1] * eigv[i];
        }
        Eigen::Matrix<double, 2, 3> AA;
        AA.col(0) = zi[1] * AMs.col(0) + zi[0] * AMs.col(2) + AMs.col(5);
        AA.col(1) = zi[1] * AMs.col(1) + zi[0] * AMs.col(4) + AMs.col(7);
        AA.col(2) = zi[1] * AMs.col(3) + zi[0] * AMs.col(6) + AMs.col(8);
        AA(0,0) = AA(0,0) - zi[2];
        AA(1,1) = AA(1,1) - zi[2];
        //AA(2,2) = AA(2,2) - zi[2];

        //Eigen::Matrix<double, 2, 1>  s = AA.leftCols(2).colPivHouseholderQr().solve(-AA.col(2));  // TODO: We only need one value.
        
        sols(0,i) = zi[0];
        //sols(1,i) = s(1);
        //std::cout << "s(1) = " << s(1) << std::endl;
        sols(1,i) = (-AA(1,2) + AA(1,0) * AA(0,2) / AA(0,0)) / (AA(1,1) - AA(1,0) * AA(0,1) / AA(0,0));
        //std::cout << "new = " << sols(1,i) << std::endl;
    }
}
} // namespace
} // namespace
