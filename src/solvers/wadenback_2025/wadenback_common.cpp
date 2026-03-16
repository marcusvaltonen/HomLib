#include <Eigen/Dense>
#include <vector>
namespace HomLib {
namespace Wadenback2025 {
double compute_residual1(
    const std::vector<Eigen::Vector2d> &x,
    double k
) {

    double b1 = x[0][1]*x[2][0] - x[0][0]*x[2][1];
    double b2 = x[1][0]*x[2][1] - x[1][1]*x[2][0];
    double b3 = x[0][0]*x[1][1] - x[0][1]*x[1][0];
    double res1 = std::abs(
        k * std::pow(x[0][0], 2) * b2 +
        k * std::pow(x[1][0], 2) * b1 +
        k * std::pow(x[0][1], 2) * b2 +
        k * std::pow(x[1][1], 2) * b1 + 
        k * std::pow(x[2][0], 2) * b3 +
        k * std::pow(x[2][1], 2) * b3 +
        b1 + b2 + b3
    );
    return res1;
}
double compute_residual2(
    const std::vector<Eigen::Vector2d> &x,
    double k
) {
    double kx00_2 = k * std::pow(x[0][0], 2);
    double kx01_2 = k * std::pow(x[0][1], 2);
    double kx10_2 = k * std::pow(x[1][0], 2);
    double kx11_2 = k * std::pow(x[1][1], 2);
    double kx20_2 = k * std::pow(x[2][0], 2);
    double kx21_2 = k * std::pow(x[2][1], 2);

    double m11 = -kx10_2*x[2][1] - kx11_2*x[2][1] + kx20_2*x[1][1] + kx21_2*x[1][1] + x[1][1] - x[2][1];
    double m12 = kx10_2*x[2][0] - kx20_2*x[1][0] - kx21_2*x[1][0] + kx11_2*x[2][0] - x[1][0] + x[2][0];
    double m13 = x[1][0]*x[2][1] - x[1][1]*x[2][0];
    double m21 = kx00_2*x[2][1] + kx01_2*x[2][1] - kx20_2*x[0][1] - kx21_2*x[0][1] - x[0][1] + x[2][1];
    double m22 = -kx00_2*x[2][0] + kx20_2*x[0][0] + kx21_2*x[0][0] - kx01_2*x[2][0] + x[0][0] - x[2][0];
    double m23 = -x[0][0]*x[2][1] + x[0][1]*x[2][0];
    double m31 = -kx00_2*x[1][1] - kx01_2*x[1][1] + k*x[0][1]*std::pow(x[1][0],2) + kx11_2*x[0][1] + x[0][1] - x[1][1];
    double m32 = kx00_2*x[1][0] - kx11_2*x[0][0] - kx11_2*x[0][0] + kx01_2*x[1][0] - x[0][0] + x[1][0];
    double m33 = x[0][0]*x[1][1] - x[0][1]*x[1][0];

    double g1 = m13 + m11*x[3][0] + m12*x[3][1];
    double g2 = m23 + m21*x[3][0] + m22*x[3][1];
    double g3 = m33 + m31*x[3][0] + m32*x[3][1];

    return g1 * g2 * g3;
}

Eigen::Matrix<double, 12, 1> compute_coeffs(
    const std::vector<Eigen::Vector2d> &x
) {
    // This writes the expression M = adj([x{1} x{2} x{3}]), such that m_ij = cij_1*k+cij_0 (note ci3_1=0, hence omitted)
    double c11_1 = -std::pow(x[1][0],2)*x[2][1] - std::pow(x[1][1],2)*x[2][1] + x[1][1]*std::pow(x[2][0],2) + x[1][1]*std::pow(x[2][1],2);
    double c11_0 = x[1][1] - x[2][1];
    double c12_1 = std::pow(x[1][0],2)*x[2][0] - x[1][0]*std::pow(x[2][0],2) - x[1][0]*std::pow(x[2][1],2) + std::pow(x[1][1],2)*x[2][0];
    double c12_0 = -x[1][0] + x[2][0];
    double c21_1 = std::pow(x[0][0],2)*x[2][1] + std::pow(x[0][1],2)*x[2][1] - x[0][1]*std::pow(x[2][0],2) - x[0][1]*std::pow(x[2][1],2);
    double c21_0 = -x[0][1] + x[2][1];
    double c22_1 = -std::pow(x[0][0],2)*x[2][0] + x[0][0]*std::pow(x[2][0],2) + x[0][0]*std::pow(x[2][1],2) - std::pow(x[0][1],2)*x[2][0];
    double c22_0 = x[0][0] - x[2][0];
    double c31_1 = -std::pow(x[0][0],2)*x[1][1] - std::pow(x[0][1],2)*x[1][1] + x[0][1]*std::pow(x[1][0],2) + x[0][1]*std::pow(x[1][1],2);
    double c31_0 = x[0][1] - x[1][1];
    double c32_1 = std::pow(x[0][0],2)*x[1][0] - x[0][0]*std::pow(x[1][0],2) - x[0][0]*std::pow(x[1][1],2) + std::pow(x[0][1],2)*x[1][0];
    double c32_0 = -x[0][0] + x[1][0];
    double c13_0 = x[1][0]*x[2][1] - x[1][1]*x[2][0];
    double c23_0 = -x[0][0]*x[2][1] + x[0][1]*x[2][0];
    double c33_0 = x[0][0]*x[1][1] - x[0][1]*x[1][0];

    // g = adj([x{1} x{2} x{3}]) * x{4}
    double r32 = x[3][0]*x[3][0] + x[3][1]*x[3][1];
    double g1_1 = c11_1*x[3][0] + c12_1*x[3][1] + c13_0*r32;
    double g1_0 = c13_0 + c11_0*x[3][0] + c12_0*x[3][1];
    double g2_1 = c21_1*x[3][0] + c22_1*x[3][1] + c23_0*r32;
    double g2_0 = c23_0 + c21_0*x[3][0] + c22_0*x[3][1];
    double g3_1 = c31_1*x[3][0] + c32_1*x[3][1] + c33_0*r32;
    double g3_0 = c33_0 + c31_0*x[3][0] + c32_0*x[3][1];

    // h = adj([x{1} x{2} x{3}]) * x{5}
    double r42 = x[4][0]*x[4][0] + x[4][1]*x[4][1];
    double h1_1 = c11_1*x[4][0] + c12_1*x[4][1] + c13_0*r42;
    double h1_0 = c13_0 + c11_0*x[4][0] + c12_0*x[4][1];
    double h2_1 = c21_1*x[4][0] + c22_1*x[4][1] + c23_0*r42;
    double h2_0 = c23_0 + c21_0*x[4][0] + c22_0*x[4][1];
    double h3_1 = c31_1*x[4][0] + c32_1*x[4][1] + c33_0*r42;
    double h3_0 = c33_0 + c31_0*x[4][0] + c32_0*x[4][1];

    // Diagonal elements of T = adj(diag(g))
    double t1_2 = g2_1*g3_1;
    double t1_1 = g2_0*g3_1 + g2_1*g3_0;
    double t1_0 = g2_0*g3_0;
    double t2_2 = g1_1*g3_1;
    double t2_1 = g1_0*g3_1 + g1_1*g3_0;
    double t2_0 = g1_0*g3_0;
    double t3_2 = g1_1*g2_1;
    double t3_1 = g1_0*g2_1 + g1_1*g2_0;
    double t3_0 = g1_0*g2_0;

    // Complete expression n = T * h, where each element of n is a monomial of degree three in the unknown distortion coefficient
    double n1_3 = h1_1*t1_2;
    double n1_2 = h1_0*t1_2 + h1_1*t1_1;
    double n1_1 = h1_0*t1_1 + h1_1*t1_0;
    double n1_0 = h1_0*t1_0;
    double n2_3 = h2_1*t2_2;
    double n2_2 = h2_0*t2_2 + h2_1*t2_1;
    double n2_1 = h2_0*t2_1 + h2_1*t2_0;
    double n2_0 = h2_0*t2_0;
    double n3_3 = h3_1*t3_2;
    double n3_2 = h3_0*t3_2 + h3_1*t3_1;
    double n3_1 = h3_0*t3_1 + h3_1*t3_0;
    double n3_0 = h3_0*t3_0;

    // return row-major order
    Eigen::Matrix<double, 12, 1> data;
    data << n1_3, n2_3, n3_3,
            n1_2, n2_2, n3_2,
            n1_1, n2_1, n3_1,
            n1_0, n2_0, n3_0;

    return data;
}


Eigen::Matrix3d homography_4pt_guo(
    const std::vector<Eigen::Vector2d> &x,
    const std::vector<Eigen::Vector2d> &y
) {
    Eigen::Matrix3d xadj;
    xadj << x[1](1)-x[2](1), x[2](0)-x[1](0), x[1](0)*x[2](1)-x[2](0)*x[1](1),
            x[2](1)-x[0](1), x[0](0)-x[2](0), x[2](0)*x[0](1)-x[0](0)*x[2](1),
            x[0](1)-x[1](1), x[1](0)-x[0](0), x[0](0)*x[1](1)-x[1](0)*x[0](1);
    Eigen::Matrix3d yadj;        
    yadj << y[1](1)-y[2](1), y[2](0)-y[1](0), y[1](0)*y[2](1)-y[2](0)*y[1](1),
            y[2](1)-y[0](1), y[0](0)-y[2](0), y[2](0)*y[0](1)-y[0](0)*y[2](1),
            y[0](1)-y[1](1), y[1](0)-y[0](0), y[0](0)*y[1](1)-y[1](0)*y[0](1);
    Eigen::Vector3d betaxl4;
    betaxl4 << xadj(0,0) * x[3](0) + xadj(0,1) * x[3](1) + xadj(0,2),
               xadj(1,0) * x[3](0) + xadj(1,1) * x[3](1) + xadj(1,2),
               xadj(2,0) * x[3](0) + xadj(2,1) * x[3](1) + xadj(2,2); 
               
    Eigen::Vector3d betaxr4;                                                                                                                          
    betaxr4 << yadj(0,0) * y[3](0) + yadj(0,1) * y[3](1) + yadj(0,2),
               yadj(1,0) * y[3](0) + yadj(1,1) * y[3](1) + yadj(1,2),
               yadj(2,0) * y[3](0) + yadj(2,1) * y[3](1) + yadj(2,2);     

    double detxl = x[1](0)*x[2](1)-x[2](0)*x[1](1)+x[2](0)*x[0](1)-x[0](0)*x[2](1)+x[0](0)*x[1](1)-x[1](0)*x[0](1);
    double t1 = detxl * betaxr4(0) / betaxl4(0);
    double t2 = detxl * betaxr4(1) / betaxl4(1);
    double t3 = detxl * betaxr4(2) / betaxl4(2);

    Eigen::Matrix3d H;
    H << t1*y[0](0), t2*y[1](0), t3*y[2](0), t1*y[0](1), t2*y[1](1), t3*y[2](1), t1, t2, t3;
    H = H * xadj;  
    return H;
}

} // namespace
} // namespace
