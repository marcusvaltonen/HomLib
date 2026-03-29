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
#include <chrono>

#include "problem_instance.hpp"
#include "radial.hpp"
#include "generate_problem_instance.hpp"

namespace HomLib {

    static const double kPI = 3.14159265358979323846;

    HomLib::ProblemInstance generate_problem_instance(const ProblemConfig &config) {
	
		HomLib::ProblemInstance instance;
		
		double fov_scale = std::tan(config.camera_fov_ / 2.0 * kPI / 180.0);

		// Random generators
		std::default_random_engine random_engine;
		random_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
		std::uniform_real_distribution<double> depth_gen(config.min_depth_, config.max_depth_);
		std::uniform_real_distribution<double> coord_gen(-fov_scale, fov_scale);
		// std::uniform_real_distribution<double> focal_gen(config.min_focal_, config.max_focal_);
		std::normal_distribution<double> direction_gen(0.0, 1.0);
		std::uniform_real_distribution<double> dist_gen(config.min_dist_, config.max_dist_);

		bool instance_generated = false;

		while (!instance_generated) {
		    Eigen::Vector3d t;
		    t.setRandom();
		    t.normalize();
		    Eigen::Matrix3d R = Eigen::Quaternion<double>::UnitRandom().toRotationMatrix();

		    // double focal_gt = focal_gen(random_engine);
		    	
	        // Point to point correspondences
	        instance.x1.clear();
	        instance.x2.clear();
	        instance.x1.reserve(config.number_points);
	        instance.x2.reserve(config.number_points);

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

					//
					Eigen::Vector2d x1h = x1.hnormalized();
		            instance.x1.push_back(x1h);
		            instance.x2.push_back(x2h);

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

		    // Distort
		    if (config.no_distortion) {
		        instance.posedata.distortion_parameter = 0.0;
		        instance.posedata.distortion_parameter2 = 0.0;
		    } else {
		        instance.posedata.distortion_parameter2 = dist_gen(random_engine);
		        if (config.one_sided) {
			    instance.posedata.distortion_parameter = 0.0;
		        } else {
			    if (config.equal) {
			        instance.posedata.distortion_parameter = instance.posedata.distortion_parameter2;
			    } else {
			        instance.posedata.distortion_parameter = dist_gen(random_engine);
			    }
		        }

		        if (!config.one_sided) { 
			    HomLib::radialdistort(instance.x1, &instance.x1, instance.posedata.distortion_parameter);
		        }
		        HomLib::radialdistort(instance.x2, &instance.x2, instance.posedata.distortion_parameter2);
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
