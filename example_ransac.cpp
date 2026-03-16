#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <Eigen/Geometry> 
#include <limits>
#include <RansacLib/ransac.h>

#include <cstdlib>
#include <ctime>
#include <time.h> 
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>
#include <chrono>
#include <sstream>
#include <iomanip>

#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_nakano_icpr_2025.hpp"
#include "get_wadenback_2025.hpp"
#include "get_kukelova_cvpr_2015.hpp"

#include "ransac_estimator.h"
#include "problem_instance.hpp"
#include "generate_problem_instance.hpp"


struct BenchmarkResults {
    std::vector<long> runtimes;
    std::vector<double> hom_err;
    std::vector<double> dist_err;   
    std::vector<std::vector<int>> inlier_history;   
};

template <typename T> void print_csv_file(std::string name, std::vector<T> v) {
    std::ofstream fd(name);
    if (fd.is_open()) {
        std::copy(v.begin(), v.end()-1, std::ostream_iterator<T>(fd, ","));
        std::copy(v.end()-1, v.end(), std::ostream_iterator<T>(fd));
        fd.close();
    } else {
        std::cout << "Unable to open file" << std::endl;
    }
}

void print_files(BenchmarkResults br, int nbr_outliers, std::string method_name) {

	std::ostringstream str;
    str << std::setw(4) << std::setfill('0') << nbr_outliers;
    print_csv_file<double>(
        std::string(method_name + "_hom_error_" + str.str() + ".csv"),
        br.hom_err
    );
    print_csv_file<double>(
        std::string(method_name + "_dist_error_" + str.str() + ".csv"),
        br.dist_err
    );
    print_csv_file<long>(
        std::string(method_name + "_timing_" + str.str() + ".csv"),
        br.runtimes
    );
    
    std::ofstream fd(std::string(method_name + "_inlier_history_" + str.str() + ".csv"));
    if (fd.is_open()) {
        for (size_t i = 0; i < br.inlier_history.size(); i++) {
            std::vector<int> v = br.inlier_history[i];
		    std::copy(v.begin(), v.end()-1, std::ostream_iterator<int>(fd, ","));
		    std::copy(v.end()-1, v.end(), std::ostream_iterator<int>(fd));
		    if (i < br.inlier_history.size() -1) {
		    	fd << std::endl;
	    	}
        }
        fd.close();
    } else {
        std::cout << "Unable to open file" << std::endl;
    }
}

template <typename Estimator> BenchmarkResults test_loransac(
    std::string name,
    Estimator* estimator,
    HomLib::ProblemConfig config,
    int nbr_outliers,
    int nbr_iter,
    int nbr_ransac_iter
) {   

    BenchmarkResults br;
    for (int k = 0; k < nbr_iter; k++) {
		HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

		std::vector<int> sample;
		int n = config.number_points;
		for( int i = 0 ; i < n ; ++i ){
		   sample.push_back(i);
		}

		auto rng = std::default_random_engine {};
        std::shuffle(std::begin(sample), std::end(sample), rng);

		for (int i = 0; i < nbr_outliers; i++) {
			Eigen::Vector2d n;
			n.setRandom();
			inst.x1[sample[i]] += n * 5000;
			n.setRandom();
			inst.x2[sample[i]] += n * 5000;
		}

		HomLib::RansacEstimator<Estimator> solver(inst.x1, inst.x2, *estimator);

		ransac_lib::LORansacOptions options;
		options.squared_inlier_threshold_ = std::pow(0.005, 2);
		options.final_least_squares_ = false;
		options.min_num_iterations_ = nbr_ransac_iter;
		options.max_num_iterations_ = nbr_ransac_iter;
		options.lo_starting_iterations_ = nbr_ransac_iter + 1;
		options.num_lsq_iterations_ = 0;
		std::srand(std::time({})); // use current time as seed for random generator
		options.random_seed_ = (unsigned int) std::rand();

		ransac_lib::LocallyOptimizedMSAC<HomLib::PoseData,
		std::vector<HomLib::PoseData>,
		HomLib::RansacEstimator<Estimator>> lomsac;
		ransac_lib::RansacStatistics ransac_stats;

		HomLib::PoseData best_model;

		auto start = std::chrono::high_resolution_clock::now();
		int num_ransac_inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);
		auto end = std::chrono::high_resolution_clock::now();
		br.runtimes.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
		
		br.hom_err.push_back(inst.hom_error(best_model.homography));
		//std::cout << "Running LOMSAC experiment with " << config.number_points
		//<< " points of which " << nbr_outliers << " are outliers for " << name << std::endl;
		//std::cout << "Homography error:   " << inst.hom_error(best_model.homography) << std::endl;
		if (config.one_sided || config.equal) {
			br.dist_err.push_back(inst.dist_error(best_model.distortion_parameter));
			//std::cout << "Dist. coeff. error: " << inst.dist_error(best_model.distortion_parameter) << std::endl;
		} else {
			br.dist_err.push_back(inst.dist_error(best_model.distortion_parameter, best_model.distortion_parameter2));
			//std::cout << "Dist. coeff. error: " << inst.dist_error(best_model.distortion_parameter, best_model.distortion_parameter2) << std::endl;
		}

		br.inlier_history.push_back(ransac_stats.inlier_history);
		/*
		std::cout << "   ... LOMSAC found " << num_ransac_inliers
		<< " inliers in " << ransac_stats.num_iterations
		<< " iterations with an inlier ratio of "
		<< ransac_stats.inlier_ratio << std::endl;
		
		std::cout << "number lo iterations " << ransac_stats.number_lo_iterations
		<< " and inlier history length " << ransac_stats.inlier_history.size() << std::endl;
		
		for (size_t j=0; j < ransac_stats.inlier_history.size(); j++) {
			std::cout << ransac_stats.inlier_history[j] << ", ";
		}

		if (ransac_stats.inlier_ratio < 0.98 * (1.0 - nbr_outliers / (double) config.number_points)) {
			std::cout << "\033[1m\033[31m FAILED!\033[0m\n" << std::endl;
		}
		*/
    }
    return br;
}


int main(int argc, char *argv[]) {
    srand((unsigned int)time(0));
    // srand(2.0);
    
    int nbr_iter = 1;
    int nbr_outliers = 0;
    double point_noise = 0.0001;
    int nbr_ransac_iter = 100;
    if (argc > 1) {
        nbr_outliers = atoi(argv[1]);
    }
    if (argc > 2) {
        nbr_iter = atoi(argv[2]);
    }
    if (argc > 3) {
        point_noise = atof(argv[3]);
    }
    if (argc > 4) {
        nbr_ransac_iter = atof(argv[4]);
    }
    std::cout << "Running tests... nbr_iter = " << nbr_iter
    	<< ", nbr_outliers = " << nbr_outliers
    	<< ", and point_noise = " << point_noise
    	<< std::endl;

    // One-sided LOMSAC tests
    HomLib::ProblemConfig config;
    config.number_points = 500;
    config.point_noise = point_noise;
    config.one_sided = true;
    config.equal = true;

    bool print_to_file = true;
    BenchmarkResults br;
    
    HomLib::FitzgibbonCVPR2001::SolverSingleSided estimator_fitzgibbon_single;
    HomLib::NakanoICPR2025::SolverSingleSided estimator_nakano_single;
    HomLib::Wadenback2025::SolverSingleSided estimator_wadenback_single;
    
    std::cout << "======== SINGLE-SIDED ========" << std::endl;

    br = test_loransac("fitzgibbon_one_sided", &estimator_fitzgibbon_single, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "fitzgibbon_one_sided");
    br = test_loransac("nakano_one_sided", &estimator_nakano_single, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "nakano_one_sided");
    br = test_loransac("wadenback_one_sided", &estimator_wadenback_single, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "wadenback_one_sided");
    
    std::cout << "======== TWO-SIDED EQUAL ========" << std::endl;
    
    config.one_sided = false;
    
    HomLib::FitzgibbonCVPR2001::SolverTwoSidedEqual estimator_fitzgibbon_two_sided_equal;
    HomLib::KukelovaCVPR2015::SolverTwoSidedEqual estimator_kukelova_two_sided_equal;
    HomLib::KukelovaCVPR2015::SolverTwoSidedEqual6Pt estimator_kukelova_two_sided_equal_6pt;
    HomLib::Wadenback2025::SolverTwoSidedEqual estimator_wadenback_two_sided_equal;
    
    br = test_loransac("fitzgibbon_equal", &estimator_fitzgibbon_two_sided_equal, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "fitzgibbon_equal");
    br = test_loransac("kukelova_two_sided_equal", &estimator_kukelova_two_sided_equal, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "kukelova_two_sided_equal");
    br = test_loransac("kukelova_two_sided_6pt_equal", &estimator_kukelova_two_sided_equal_6pt, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "kukelova_two_sided_6pt_equal");
    br = test_loransac("wadenback_two_sided_equal", &estimator_wadenback_two_sided_equal, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "wadenback_two_sided_equal");
    
    std::cout << "======== TWO-SIDED ========" << std::endl;
    
    config.equal = false;
    
    HomLib::KukelovaCVPR2015::SolverTwoSided estimator_kukelova_two_sided;
    HomLib::KukelovaCVPR2015::SolverTwoSided6Pt estimator_kukelova_two_sided_6pt;
    HomLib::Wadenback2025::SolverTwoSided estimator_wadenback_two_sided;

    br = test_loransac("kukelova_two_sided", &estimator_kukelova_two_sided, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "kukelova_two_sided");
    br = test_loransac("kukelova_two_sided_6pt", &estimator_kukelova_two_sided_6pt, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "kukelova_two_sided_6pt");
    br = test_loransac("wadenback_two_sided", &estimator_wadenback_two_sided, config, nbr_outliers, nbr_iter, nbr_ransac_iter);
    if (print_to_file)
        print_files(br, nbr_outliers, "wadenback_two_sided");
    
    // Check normalization
    /*
    std::cout << "======== Test... no normalization of image coordinates ========" << std::endl;
    estimator_wadenback_two_sided.normalize_image_coord = false;
    test_loransac("Wadenback two-sided", &estimator_wadenback_two_sided, config);
    */
    
    
    
}
