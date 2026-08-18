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
#include <random>
#include <cmath>
#include <chrono>

#include "problem_instance.hpp"
#include "radial.hpp"
#include "generate_problem_instance.hpp"
#include "affine2sift.hpp"

namespace HomLib {

    static const double kPI = 3.14159265358979323846;
    
    static Eigen::Matrix2d affineFromHomography( const Eigen::Matrix3d &H, const Eigen::Vector2d &x, const Eigen::Vector2d &y , double k)
    {
        // x are the distorted coeffs
        double h1 = H(0,0), h2 = H(0,1), h3 = H(0,2),
               h4 = H(1,0), h5 = H(1,1), h6 = H(1,2),
               h7 = H(2,0), h8 = H(2,1), h9 = H(2,2);
        double u1 = x(0), v1 = x(1),
               u2 = y(0), v2 = y(1);
        double dist_fact = k*(u1*u1 + v1*v1) + 1;
        double s = h7*u1 + h8*v1 + h9*(dist_fact);
         // Note that
        // u2 = (h1*u1 + h2*v1 + h3*(dist_fact))/(s)
        // v2 = (h4*u1 + h5*v1 + h6*(dist_fact))/(s)
        Eigen::Matrix2d A;
        A << (h1 + 2*h3*k*u1)/(s) - ((h7 + 2*h9*k*u1)*(h1*u1 + h2*v1 + h3*(dist_fact)))/(s*s),
             (h2 + 2*h3*k*v1)/(s) - ((h8 + 2*h9*k*v1)*(h1*u1 + h2*v1 + h3*(dist_fact)))/(s*s),
             (h4 + 2*h6*k*u1)/(s) - ((h7 + 2*h9*k*u1)*(h4*u1 + h5*v1 + h6*(dist_fact)))/(s*s),
             (h5 + 2*h6*k*v1)/(s) - ((h8 + 2*h9*k*v1)*(h4*u1 + h5*v1 + h6*(dist_fact)))/(s*s);
        

        return A;
    }

    HomLib::ProblemInstance generate_problem_instance(const ProblemConfig &config) {
	
		HomLib::ProblemInstance instance;
		
		double fov_scale = std::tan(config.camera_fov_ / 2.0 * kPI / 180.0);

		// Random generators
		std::default_random_engine random_engine;
		random_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
		std::uniform_real_distribution<double> depth_gen(config.min_depth_, config.max_depth_);
		std::uniform_real_distribution<double> coord_gen(-fov_scale, fov_scale);
		std::uniform_real_distribution<double> focal_gen(config.min_focal_, config.max_focal_);
		std::normal_distribution<double> direction_gen(0.0, 1.0);
		std::uniform_real_distribution<double> dist_gen(config.min_dist_, config.max_dist_);

		bool instance_generated = false;

		while (!instance_generated) {
		    Eigen::Vector3d t;
		    t.setRandom();
		    t.normalize();
		    Eigen::Matrix3d R = Eigen::Quaternion<double>::UnitRandom().toRotationMatrix();

		    double focal_gt = focal_gen(random_engine);
		    	
	        // Point to point correspondences
	        instance.x1.clear();
	        instance.x2.clear();
	        instance.A.clear();
	        instance.ori.clear();
	        instance.x1.reserve(config.number_points);
	        instance.x2.reserve(config.number_points);
	        instance.A.reserve(config.number_points);
	        instance.ori.reserve(config.number_points);

		    // Generate plane
		    Eigen::Vector3d n;
		    n << direction_gen(random_engine), direction_gen(random_engine), direction_gen(random_engine);
		    n.normalize();

		    // Choose depth of plane such that center point of image 1 is at depth d
		    double d_center = depth_gen(random_engine);
		    double alpha = d_center / n(2);
		    // plane is n'*X = alpha

		    // ground truth homography
		    instance.posedata.homography = alpha * R + t * n.transpose();
		    
		    // Distort
		    switch (config.distortion) {
        	    case HomLib::DistortionCase::NO_DISTORTION:
                    instance.posedata.distortion_parameter = 0.0;
		            instance.posedata.distortion_parameter2 = 0.0;
		            break;
        	    case HomLib::DistortionCase::ONE_SIDED_LEFT:
                    instance.posedata.distortion_parameter = 0.0;
		            instance.posedata.distortion_parameter2 = dist_gen(random_engine);
		            break;
        	    case HomLib::DistortionCase::ONE_SIDED_RIGHT:
                    instance.posedata.distortion_parameter = dist_gen(random_engine);
		            instance.posedata.distortion_parameter2 = 0.0;
		            break;
        	    case HomLib::DistortionCase::TWO_SIDED_EQUAL:
                    instance.posedata.distortion_parameter = dist_gen(random_engine);
		            instance.posedata.distortion_parameter2 = instance.posedata.distortion_parameter;
		            break;
        	    case HomLib::DistortionCase::TWO_SIDED: 
        	        instance.posedata.distortion_parameter = dist_gen(random_engine);
		            instance.posedata.distortion_parameter2 = dist_gen(random_engine);
		            break;
    	    }

		    bool failed_instance = false;
		    for (int j = 0; j < config.number_points; ++j) {
		        bool point_okay = false;
		        for (int trials = 0; trials < 10; ++trials) {
		            Eigen::Vector3d x1{coord_gen(random_engine), coord_gen(random_engine), 1.0};
		            x1.normalize();
		            Eigen::Vector3d X;

		            // compute depth
		            double lambda = alpha / n.dot(x1);
		            X = x1 * lambda;
		            // Map into second image
		            X = R * X + t;

		            Eigen::Vector3d x2 = X.normalized();

		            // Check cheirality
		            if (x2(2) < 0 || lambda < 0) {
		                // try to generate another point
		                continue;
		            }

		            // Check FoV of second camera
		            Eigen::Vector2d x2h = x2.hnormalized();
		            if (x2h(0) < -fov_scale || x2h(0) > fov_scale || x2h(1) < -fov_scale || x2h(1) > fov_scale) {
		                // try to generate another point
		                continue;
		            }
			        Eigen::Vector2d x1h = x1.hnormalized();
			        
			        // Distort
			        Eigen::Vector2d x1hd, x2hd;
			        x1hd = HomLib::radialdistort(x1h, instance.posedata.distortion_parameter);
		            x2hd = HomLib::radialdistort(x2h, instance.posedata.distortion_parameter2);
		            // calculate affine from homography
		            // This assumes DistortionCase.NO_DISTORTION or DistortionCase.ONE_SIDED_RIGHT
		            Eigen::Matrix2d A = affineFromHomography(instance.posedata.homography, x1hd, x2hd, instance.posedata.distortion_parameter);
		            
		            // check if determinant is positive
		            if ( A.determinant() < 0 ) continue;

					//
					double s_1, c_1, s_2, c_2, q;
					affine2sift(A, s_1, c_1, s_2, c_2, q);
                    /*  ONLY WORKS FOR RIGHT-SIDED AND NO_DISTORTION
                    double h_1 = instance.posedata.homography(0,0),
                        h_2 = instance.posedata.homography(0,1),
                        h_3 = instance.posedata.homography(0,2),
                        h_4 = instance.posedata.homography(1,0),
                        h_5 = instance.posedata.homography(1,1),
                        h_6 = instance.posedata.homography(1,2),
                        h_7 = instance.posedata.homography(2,0),
                        h_8 = instance.posedata.homography(2,1),
                        h_9 = instance.posedata.homography(2,2);
                    double u_1 = x1hd[0],
                        v_1 = x1hd[1],
                        u_2 = x2hd[0],
                        v_2 = x2hd[1];
                    lambda = instance.posedata.distortion_parameter;

					double res = -h_1*s_2*c_1 - h_2*s_1*s_2 + h_4*c_1*c_2 + h_5*s_1*c_2 + h_7*u_2*s_2*c_1 - h_7*v_2*c_1*c_2 + h_8*u_2*s_1*s_2
                        -h_8*v_2*s_1*c_2 -2*h_3*u_1*s_2*c_1*lambda - 2*h_3*v_1*s_1*s_2*lambda + 2*h_6*u_1*c_1*c_2*lambda + 2*h_6*v_1*s_1*c_2*lambda
                        + 2*h_9*u_1*u_2*s_2*c_1*lambda - 2*h_9*u_1*v_2*c_1*c_2*lambda + 2*h_9*v_1*u_2*s_1*s_2*lambda - 2*h_9*v_1*v_2*s_1*c_2*lambda;
                    */

					Eigen::Vector2d ori;
					ori[0] = std::atan2(s_1, c_1);
					ori[1] = std::atan2(s_2, c_2);

		            instance.x1.push_back(x1hd);
		            instance.x2.push_back(x2hd);
		            instance.A.push_back(A);
		            instance.ori.push_back(ori);

		            point_okay = true;
		            break;
		        }
		        if (!point_okay) {
		            failed_instance = true;
		            break;
		        }
		    }
		    if (failed_instance) {
		        continue;
		    }
		    
		    // Focal length		
			//instance.x1 *= focal_gt;
			//instance.x2 *= focal_gt;
			//instance.posedata  // mod homography and params
			

			// Add noise
		    if (config.point_noise > 0) {
		        std::normal_distribution<double> normal;
		        for (int i = 0; i < config.number_points; i++) {
		            instance.x1[i](0) += normal(random_engine) * config.point_noise;
		            instance.x1[i](1) += normal(random_engine) * config.point_noise;
		            instance.x2[i](0) += normal(random_engine) * config.point_noise;
		            instance.x2[i](1) += normal(random_engine) * config.point_noise;
		        }
		    }
		    instance_generated = true;
		}
		return instance;
	}
}  // namespace HomLib
