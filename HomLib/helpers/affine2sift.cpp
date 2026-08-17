#include <vector>
#include <numeric>
#include <Eigen/Dense>





using namespace Eigen;



void affine2sift_solver(Eigen::MatrixXd const& A, double q1, std::vector<double>* r1, std::vector<double>* r2)
{
    // Compute coefficients
double const* A_data = A.data();
double const A11 = A_data[0];
double const A21 = A_data[1];
double const A12 = A_data[2];
double const A22 = A_data[3];
    VectorXd coeffs(14);
  double _t2_ = A12*2.0;
  double _t3_ = A22*2.0;
  double _t4_ = q1*2.0;
  double _t5_ = -A11;
  double _t6_ = -A21;
  double _t7_ = -q1;
  double _t8_ = -_t4_;
  coeffs[0] = q1+_t5_;
  coeffs[1] = _t5_+_t7_;
  coeffs[2] = _t2_;
  coeffs[3] = _t2_;
  coeffs[4] = A11+q1;
  coeffs[5] = A11+_t7_;
  coeffs[6] = _t6_;
  coeffs[7] = _t8_;
  coeffs[8] = _t6_;
  coeffs[9] = _t3_;
  coeffs[10] = _t3_;
  coeffs[11] = A21;
  coeffs[12] = _t8_;
  coeffs[13] = A21;


    // Setup elimination template
    static const int coeffs0_ind[] = { 0,6,0,6,7,0,2,6,9,1,7,8,1,8,1,3,7,8,10,0,2,6,7,9,2,4,9,11 };
    static const int coeffs1_ind[] = { 5,13,3,5,10,13,1,3,8,10,3,5,10,12,13,2,4,9,11,12,5,12,13,4,11,12,4,11 };


static const int C0_ind[] = {3,7,10,14,15,17,19,21,23,27,30,31,34,38,41,43,44,45,47,48,50,52,53,54,57,59,61,63};

static const int C1_ind[] = {0,4,8,10,12,14,16,18,20,22,25,27,29,30,31,32,34,36,38,39,41,44,45,48,52,53,57,61};

MatrixXd C0 = MatrixXd::Zero(8,8);
MatrixXd C1 = MatrixXd::Zero(8,8);
for (int i = 0; i < 28; i++) {
    C0(C0_ind[i]) = coeffs(coeffs0_ind[i]);
}

for (int i = 0; i < 28; i++) {
    C1(C1_ind[i]) = coeffs(coeffs1_ind[i]);
}

/*
Eigen::Matrix<double,8,4> b = Eigen::Matrix<double,8,4>::Zero();
b(4,0) = -1;
b(5,1) = -1;
b(6,2) = -1;
b(7,3) = -1;
Eigen::Matrix<double,8,4> alpha = C0.transpose().fullPivLu().solve(b);
Eigen::Matrix<double,12,4> RR;
RR << alpha.transpose()*C1, Eigen::Matrix<double,8,8>::Identity();
//AM_ind = [6,7,1,2,3,8,9,4];
//AM = RR(AM_ind,:);
Eigen::Matrix<double,8,8> AM;
AM << RR.col(5), RR.col(6), RR.col(0), RR.col(1), RR.col(2), RR.col(7), RR.col(8), RR.col(3);
Eigen::EigenSolver< Eigen::Matrix<double,8,8> > AMsolver(AM);
Eigen::MatrixXcd V = AMsolver.eigenvectors();
V = V.array() * (Eigen::Matrix<double,8,1>::Ones()*V.col(0)).array();
Eigen::VectorXcd r1c = AMsolver.eigenvalues();
Eigen::VectorXcd r2c = V.row(5);
*/

//[V,D] = eig(AM);
//V = V ./ (ones(size(V,1),1)*V(1,:));
//sols(1,:) = diag(D).';
//sols(2,:) = V(6,:);

MatrixXd C12 = C0.fullPivLu().solve(C1);



    // Setup action matrix
    Matrix<double,12, 8> RR;
    RR << -C12.bottomRows(4), Matrix<double,8,8>::Identity(8, 8);

    static const int AM_ind[] = { 5,6,0,1,2,7,8,3 };
    Matrix<double, 8, 8> AM;
    for (int i = 0; i < 8; i++) {
        AM.row(i) = RR.row(AM_ind[i]);
    }

    MatrixXcd sols(2, 8);
    sols.setZero();

    // Solve eigenvalue problem
    EigenSolver<Matrix<double, 8, 8> > es(AM);
    ArrayXcd D = es.eigenvalues();
    ArrayXXcd V = es.eigenvectors();

    V = (V / V.row(0).array().replicate(8, 1)).eval();


        sols.row(0) = D.transpose().array();
    sols.row(1) = V.row(5).array();




    Eigen::VectorXcd r1c = sols.row(0);
    Eigen::VectorXcd r2c = sols.row(1);
    int nsols = r1c.size();
    for (int isol = 0; isol < nsols; ++isol) {
        if ( r1c(isol).imag() == 0 && r2c(isol).imag() == 0 )
        {
            r1->push_back(r1c(isol).real());
            r2->push_back(r2c(isol).real());
        }
    }
}

// Action =
// Quotient ring basis (V) =  r1^2, r1*r2, r1*r2^2, r2, r2^2, r2^3
// Available monomials (RR*V) = r1^2*r2, r1^2*r2^2, r1*r2^3, 1, r1, r1^2, r1*r2, r1*r2^2, r2, r2^2, r2^3

void affine2sift(const Eigen::Matrix2d &A, double &s1, double &c1, double &s2, double &c2, double &q )
{
    q = sqrt(A.determinant());
    std::vector<double> r1solns, r2solns;
    affine2sift_solver(A, q, &r1solns, &r2solns);

    double r1 = r1solns[0];
    double r2 = r2solns[0];
    c1 = (1-r1*r1)/(1+r1*r1);
    s1 = (2*r1)/(1+r1*r1);
    c2 = (1-r2*r2)/(1+r2*r2);
    s2 = (2*r2)/(1+r2*r2);

    // check residuals
    //double res1 = c1*s2*A(0,0) + s1*s2*A(0,1) - c1*c2*A(1,0) - c2*s1*A(1,1);
    //double res2 = A(0,1)*A(1,0)-A(0,0)*A(1,1)+q*q;
    //double res3 = A(0,0)*c1 + A(0,1)*s1 - c2*q;
    //double res4 = A(1,0)*c1 + A(1,1)*s1 - s2*q;
}
