#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "posedata.hpp"

static const double SMALL_NUMBER = 1e-8;
static const double TOL_CONVERGENCE = 1e-10;
static const double INITIAL_LM_DAMP = 1e-4;
static const int MAX_ITER = 10;


namespace HomLib {
// Non-linear refinement of transfer error |x2 - pi(H*x1)|^2 with no distortion (k = 0), parameterized by fixing H(2,2) = 1
void refinement_unsided(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    int num_points = static_cast<int>(x1.size());
    int num_residuals = 2 * num_points;
    int num_parameters = 8; // Only homography parameters
    double lm_damping = INITIAL_LM_DAMP;

    // Order for jacobian is: h11, h21, h31, h12, h22, h32, h13, h23
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J(num_residuals, num_parameters);
    J.setZero();
    Eigen::Matrix<double, Eigen::Dynamic, 1> residuals(num_residuals, 1);
    Eigen::Matrix<double, Eigen::Dynamic, 1> delta(num_parameters, 1);
    residuals.setZero();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian;
    Eigen::Matrix<double, Eigen::Dynamic, 1> gradient;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        Eigen::Matrix3d H = p.homography;
        H /= H(2, 2);

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);

        for (size_t i = 0; i < x1.size(); i++) {
            const double x1_0 = x1[i](0);
            const double x1_1 = x1[i](1);
            const double x2_0 = x2[i](0);
            const double x2_1 = x2[i](1);

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;

            // Residuals (no distortion)
            residuals(2 * i + 0) = z0 - x2_0;
            residuals(2 * i + 1) = z1 - x2_1;

            // Jacobian w.r.t. homography elements
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;               // h11
            // J(2 * i + 0, 1) = 0.0;                          // h21
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;         // h31
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;               // h12
            // J(2 * i + 0, 4) = 0.0;                          // h22
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;         // h32
            J(2 * i + 0, 6) = 1.0 * inv_Hx1_2;                // h13
            // J(2 * i + 0, 7) = 0.0;                          // h23

            // Second row
            // J(2 * i + 1, 0) = 0.0;                          // h11
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;               // h21
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;         // h31
            // J(2 * i + 1, 3) = 0.0;                          // h12
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;               // h22
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;         // h32
            // J(2 * i + 1, 6) = 0.0;                          // h13
            J(2 * i + 1, 7) = 1.0 * inv_Hx1_2;                // h23
        }

        if (residuals.norm() < TOL_CONVERGENCE)
            break;

        hessian = J.transpose() * J;
        hessian.diagonal().array() += lm_damping; // LM
        gradient = -J.transpose() * residuals;

        if (gradient.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
            break;

        delta = hessian.ldlt().solve(gradient);

        Eigen::Matrix3d H_new = H;
        Eigen::Map<Eigen::Matrix<double, 8, 1>>(H_new.data()) += delta.head(8);
        p.homography = H_new;

        if (delta.array().abs().maxCoeff() < SMALL_NUMBER)
            break;
        lm_damping = std::max(1e-8, lm_damping / 10.0);
    }
}
// Non-linear refinement of transfer error |x2(k) - pi(H*x1)|^2, parameterized by fixing H(2,2) = 1
void refinement_onesided(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    int n_pts = x1.size();
	int n_res = 2 * n_pts;
	int n_params = 9;
	double lm_damp = INITIAL_LM_DAMP;
	
	// Order for jacobian is: h11, h21, h31, h12, h22, h32, h13, h23, k
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J(n_res, n_params);
	J.setZero();
	Eigen::Matrix<double, Eigen::Dynamic, 1> res(n_res, 1);
	Eigen::Matrix<double, Eigen::Dynamic, 1> dx(n_params, 1);
	res.setZero();
	
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hess;
	Eigen::Matrix<double, Eigen::Dynamic, 1> g;
	
    //std::cout << "allocation done\n";
	
	for (int iter = 0; iter < MAX_ITER; iter++) {
	    Eigen::Matrix3d H = p.homography;
	    H /= H(2,2);
        const double k = p.distortion_parameter;

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);
        
        //std::cout << "entering first row of iter"<<iter <<"\n";

        for (size_t i = 0; i < x1.size(); i++) {
            const double x1_0 = x1[i](0);
            const double x1_1 = x1[i](1);
            const double x2_0 = x2[i](0);
            const double x2_1 = x2[i](1);

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;
            
            // Residual
            const double r22 = (x2_0 * x2_0 + x2_1 * x2_1);
            const double dist_factor = 1.0 / (1.0 + k * r22);
            res(2 * i + 0) = z0 - x2_0 * dist_factor;
            res(2 * i + 1) = z1 - x2_1 * dist_factor;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = 1.0 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;
            J(2 * i + 0, 8) = x2_0 * r22 * dist_factor * dist_factor;  // Dist coeff.
            
            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = 1.0 * inv_Hx1_2;
            J(2 * i + 1, 8) = x2_1 * r22 * dist_factor * dist_factor;  // Dist coeff.
            
            //std::cout << "iter=" << iter << "row=" << i << " J=\n" << J << "\n";
        }
        
        //std::cout << "res=\n" << res << "\n";
        //std::cout << "J=\n" << J << "\n";

	    if (res.norm() < TOL_CONVERGENCE)
		    break;

	    Hess = J.transpose()*J;
	    Hess.diagonal().array() += lm_damp; // LM dampening
	    g = -J.transpose()*res;


	    if (g.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
		    break;

	    //std::cout << "iter=" << iter << " res=" << res.squaredNorm() << ", g="<< g.squaredNorm() << "\n";
	    //std::cout << res << "\n";

	    dx = Hess.ldlt().solve(g);
	    
	    Eigen::Matrix3d H_new = H;
	    Eigen::Map<Eigen::Matrix<double, 8, 1>>(H_new.data()) += dx.head(8);
        p.homography = H_new;
	    p.distortion_parameter += dx(8);

	    if (dx.array().abs().maxCoeff() < SMALL_NUMBER)
		    break;
	    lm_damp = std::max(1e-8, lm_damp / 10.0);
	}
}
// Non-linear refinement of transfer error |x2(k) - pi(H*x1(k))|^2, parameterized by fixing H(2,2) = 1
void refinement_twosided_equal(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    int n_pts = x1.size();
	int n_res = 2 * n_pts;
	int n_params = 9;
	double lm_damp = INITIAL_LM_DAMP;
	
	// Order for jacobian is: h11, h21, h31, h12, h22, h32, h13, h23, k
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J(n_res, n_params);
	J.setZero();
	Eigen::Matrix<double, Eigen::Dynamic, 1> res(n_res, 1);
	Eigen::Matrix<double, Eigen::Dynamic, 1> dx(n_params, 1);
	res.setZero();
	
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hess;
	Eigen::Matrix<double, Eigen::Dynamic, 1> g;
	
    //std::cout << "allocation done\n";
	
	for (int iter = 0; iter < MAX_ITER; iter++) {
	    Eigen::Matrix3d H = p.homography;
	    H /= H(2,2);
        const double k = p.distortion_parameter;

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);
        
        //std::cout << "entering first row of iter"<<iter <<"\n";

        for (size_t i = 0; i < x1.size(); i++) {          
            double x1_0 = x1[i](0);
            double x1_1 = x1[i](1);
            double x2_0 = x2[i](0);
            double x2_1 = x2[i](1);
            
            const double r12 = (x1_0 * x1_0 + x1_1 * x1_1);
            const double inv_dist_factor1 = 1.0 + k * r12;

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2 * inv_dist_factor1;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2 * inv_dist_factor1;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2 * inv_dist_factor1);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;
            
            // Residual
            const double r22 = (x2_0 * x2_0 + x2_1 * x2_1);
            const double dist_factor2 = 1.0 / (1.0 + k * r22);
            res(2 * i + 0) = z0 - x2_0 * dist_factor2;
            res(2 * i + 1) = z1 - x2_1 * dist_factor2;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = 1.0 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;
            J(2 * i + 0, 8) = x2_0 * r22 * dist_factor2 * dist_factor2 + H0_2 * r12 * inv_Hx1_2 - z0 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff.
            
            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = 1.0 * inv_Hx1_2;
            J(2 * i + 1, 8) = x2_1 * r22 * dist_factor2 * dist_factor2 + H1_2 * r12 * inv_Hx1_2 - z1 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff.
        }
        
        // std::cout << "res=\n" << res << "\n";
        // std::cout << "J=\n" << J << "\n";

	    if (res.norm() < TOL_CONVERGENCE)
		    break;

	    Hess = J.transpose()*J;
	    Hess.diagonal().array() += lm_damp; // LM dampening
	    g = -J.transpose()*res;


	    if (g.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
		    break;

	    //std::cout << "iter=" << iter << " res=" << res.squaredNorm() << ", g="<< g.squaredNorm() << "\n";
	    //std::cout << res << "\n";

	    dx = Hess.ldlt().solve(g);
	    
	    Eigen::Matrix3d H_new = H;
	    Eigen::Map<Eigen::Matrix<double, 8, 1>>(H_new.data()) += dx.head(8);
        p.homography = H_new;
	    p.distortion_parameter += dx(8);

	    if (dx.array().abs().maxCoeff() < SMALL_NUMBER)
		    break;
	    lm_damp = std::max(1e-8, lm_damp / 10.0);
	}
}
// Non-linear refinement of transfer error |x2(k2) - pi(H*x1(k1))|^2, parameterized by fixing H(2,2) = 1
void refinement_twosided(
    const std::vector<Eigen::Vector2d> &x1,
    const std::vector<Eigen::Vector2d> &x2,
    HomLib::PoseData &p
) {
    int n_pts = x1.size();
	int n_res = 2 * n_pts;
	int n_params = 10;
	double lm_damp = INITIAL_LM_DAMP;
	
	// Order for jacobian is: h11, h21, h31, h12, h22, h32, h13, h23, k
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J(n_res, n_params);
	J.setZero();
	Eigen::Matrix<double, Eigen::Dynamic, 1> res(n_res, 1);
	Eigen::Matrix<double, Eigen::Dynamic, 1> dx(n_params, 1);
	res.setZero();
	
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hess;
	Eigen::Matrix<double, Eigen::Dynamic, 1> g;
	
    //std::cout << "allocation done\n";
	
	for (int iter = 0; iter < MAX_ITER; iter++) {
	    Eigen::Matrix3d H = p.homography;
	    H /= H(2,2);
        const double k1 = p.distortion_parameter;
        const double k2 = p.distortion_parameter2;

        const double H0_0 = H(0, 0), H0_1 = H(0, 1), H0_2 = H(0, 2);
        const double H1_0 = H(1, 0), H1_1 = H(1, 1), H1_2 = H(1, 2);
        const double H2_0 = H(2, 0), H2_1 = H(2, 1), H2_2 = H(2, 2);
        
        //std::cout << "entering first row of iter"<<iter <<"\n";

        for (size_t i = 0; i < x1.size(); i++) {          
            double x1_0 = x1[i](0);
            double x1_1 = x1[i](1);
            double x2_0 = x2[i](0);
            double x2_1 = x2[i](1);
            
            const double r12 = (x1_0 * x1_0 + x1_1 * x1_1);
            const double inv_dist_factor1 = 1.0 + k1 * r12;

            const double Hx1_0 = H0_0 * x1_0 + H0_1 * x1_1 + H0_2 * inv_dist_factor1;
            const double Hx1_1 = H1_0 * x1_0 + H1_1 * x1_1 + H1_2 * inv_dist_factor1;
            const double inv_Hx1_2 = 1.0 / (H2_0 * x1_0 + H2_1 * x1_1 + H2_2 * inv_dist_factor1);

            const double z0 = Hx1_0 * inv_Hx1_2;
            const double z1 = Hx1_1 * inv_Hx1_2;
            
            // Residual
            const double r22 = (x2_0 * x2_0 + x2_1 * x2_1);
            const double dist_factor2 = 1.0 / (1.0 + k2 * r22);
            res(2 * i + 0) = z0 - x2_0 * dist_factor2;
            res(2 * i + 1) = z1 - x2_1 * dist_factor2;

            // Jacobian w.r.t. first element
            J(2 * i + 0, 0) = x1_0 * inv_Hx1_2;
            // J(2 * i + 0, 1) = 0.0;
            J(2 * i + 0, 2) = -x1_0 * z0 * inv_Hx1_2;
            J(2 * i + 0, 3) = x1_1 * inv_Hx1_2;
            // J(2 * i + 0, 4) = 0.0;
            J(2 * i + 0, 5) = -x1_1 * z0 * inv_Hx1_2;
            J(2 * i + 0, 6) = 1.0 * inv_Hx1_2;
            // J(2 * i + 0, 7) = 0.0;
            J(2 * i + 0, 8) = H0_2 * r12 * inv_Hx1_2 - z0 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff. k1
            J(2 * i + 0, 9) = x2_0 * r22 * dist_factor2 * dist_factor2;  // Dist coeff. k2
            
            // Jacobian w.r.t. second element
            // J(2 * i + 1, 0) = 0.0;
            J(2 * i + 1, 1) = x1_0 * inv_Hx1_2;
            J(2 * i + 1, 2) = -x1_0 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 3) = 0.0;
            J(2 * i + 1, 4) = x1_1 * inv_Hx1_2;
            J(2 * i + 1, 5) = -x1_1 * z1 * inv_Hx1_2;
            // J(2 * i + 1, 6) = 0.0;
            J(2 * i + 1, 7) = 1.0 * inv_Hx1_2;
            J(2 * i + 1, 8) = H1_2 * r12 * inv_Hx1_2 - z1 * H2_2 * r12 * inv_Hx1_2;  // Dist coeff. k1
            J(2 * i + 1, 9) = x2_1 * r22 * dist_factor2 * dist_factor2;  // Dist coeff. k2
        }
        
        // std::cout << "res=\n" << res << "\n";
        // std::cout << "J=\n" << J << "\n";

	    if (res.norm() < TOL_CONVERGENCE)
		    break;

	    Hess = J.transpose()*J;
	    Hess.diagonal().array() += lm_damp; // LM dampening
	    g = -J.transpose()*res;


	    if (g.cwiseAbs().maxCoeff() < TOL_CONVERGENCE)
		    break;

	    //std::cout << "iter=" << iter << " res=" << res.squaredNorm() << ", g="<< g.squaredNorm() << "\n";
	    //std::cout << res << "\n";

	    dx = Hess.ldlt().solve(g);
	    
	    Eigen::Matrix3d H_new = H;
	    Eigen::Map<Eigen::Matrix<double, 8, 1>>(H_new.data()) += dx.head(8);
        p.homography = H_new;
	    p.distortion_parameter += dx(8);
	    p.distortion_parameter2 += dx(9);

	    if (dx.array().abs().maxCoeff() < SMALL_NUMBER)
		    break;
	    lm_damp = std::max(1e-8, lm_damp / 10.0);
	}
}
}

