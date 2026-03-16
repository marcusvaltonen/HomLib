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
#include <vector>
#include <numeric>
#include <chrono>  // NOLINT [build/c++11]
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_kukelova_cvpr_2015.hpp"
#include "get_nakano_icpr_2025.hpp"
// #include "get_valtonenornhag_icpr_2020.hpp"
// #include "get_valtonenornhag_wacv_2021.hpp"
#include "get_wadenback_3dv_2026.hpp"
#include "problem_instance.hpp"
#include "generate_problem_instance.hpp"
#include "posedata.hpp"
#include "refinement.hpp"
#include "radial.hpp"

#include "wadenback_common.hpp"
#include "PoseLib/solvers/homography_4pt.h"

struct BenchmarkResults {
    std::vector<long> runtimes;
    std::vector<double> hom_err;
    std::vector<double> dist_err;   
};

struct SolverKukelova {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::KukelovaCVPR2015::get(inst.x1, inst.x2, false);
    }   
};
struct SolverKukelovaEqual {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::KukelovaCVPR2015::get(inst.x1, inst.x2, true);
    }
};
struct SolverKukelova6pt {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::KukelovaCVPR2015::get_6pt(inst.x1, inst.x2, false);
    }   
};
struct SolverKukelova6ptEqual {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::KukelovaCVPR2015::get_6pt(inst.x1, inst.x2, true);
    }   
};
struct SolverFitzgibbon {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::FitzgibbonCVPR2001::get(inst.x1, inst.x2);
    }   
};
struct SolverFitzgibbonSingle {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::FitzgibbonCVPR2001::get_single(inst.x1, inst.x2);
    }   
};
struct SolverNakano {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::NakanoICPR2025::get(inst.x1, inst.x2, false);
    }   
};
struct SolverWadenbackDouble {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::Wadenback2025::get_double_sided(inst.x1, inst.x2, false);
    }   
};
struct SolverWadenbackDoubleEqual {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::Wadenback2025::get_double_sided_equal(inst.x1, inst.x2, false);
    }   
};
struct SolverWadenbackOne {
    static inline std::vector<HomLib::PoseData> solve(const HomLib::ProblemInstance inst) {
        return HomLib::Wadenback2025::get_one_sided(inst.x1, inst.x2, false);
    }   
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

void print_files(BenchmarkResults br, double point_noise, std::string method_name) {
    print_csv_file<double>(
        std::string(method_name + "_hom_error_" + std::to_string(point_noise) + ".csv"),
        br.hom_err
    );
    print_csv_file<double>(
        std::string(method_name + "_dist_error_" + std::to_string(point_noise) + ".csv"),
        br.dist_err
    );
    print_csv_file<long>(
        std::string(method_name + "_timing_" + std::to_string(point_noise) + ".csv"),
        br.runtimes
    );
}

template <typename Solver> BenchmarkResults benchmark_solver(HomLib::ProblemConfig config, int nbr_iter) {
    BenchmarkResults br;
    std::cout << "Not found: ";
    for (int i = 0; i < nbr_iter; i++) {
        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<HomLib::PoseData> pd = Solver::solve(inst);
        auto end = std::chrono::high_resolution_clock::now();
        br.runtimes.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
        
        // Check solutions
        std::vector<double> hom_err_put;
        std::vector<double> dist_err_put;
        for (size_t j = 0; j < pd.size(); j++)
        {
            hom_err_put.push_back(inst.hom_error(pd[j].homography));
            if (config.one_sided || config.equal) {
                dist_err_put.push_back(inst.dist_error(pd[j].distortion_parameter));
            } else {
                dist_err_put.push_back(inst.dist_error(pd[j].distortion_parameter, pd[j].distortion_parameter2));
            }
        }
        
        // Handle case when no solution is found
        if ((hom_err_put.size() > 0) && (dist_err_put.size() > 0)) {
            br.hom_err.push_back(std::log10(*std::min_element(hom_err_put.begin(),hom_err_put.end())));
            br.dist_err.push_back(std::log10(*std::min_element(dist_err_put.begin(),dist_err_put.end())));
        } else {
            //br.hom_err.push_back(0.0);
            //br.dist_err.push_back(0.0);
            std::cout << "x";
        }
    }
    std::cout << "\n";
    std::sort(br.runtimes.begin(), br.runtimes.end());
    long runtime_ns = br.runtimes[br.runtimes.size() / 2];
    
    std::sort(br.hom_err.begin(), br.hom_err.end());
    double hom_err_median = br.hom_err[br.hom_err.size() / 2];
    double hom_err_mean = std::accumulate(br.hom_err.begin(), br.hom_err.end(), 0.0) / br.hom_err.size();
    
    std::sort(br.dist_err.begin(), br.dist_err.end());
    double dist_err_median = br.dist_err[br.dist_err.size() / 2];
    double dist_err_mean = std::accumulate(br.dist_err.begin(), br.dist_err.end(), 0.0) / br.dist_err.size();
    
    std::cout << "\tMedian execution time: " << runtime_ns << " ns" << std::endl;
    std::cout << "\tHomography error: " << hom_err_median << " (median), " << hom_err_mean << " (mean)" << std::endl;
    std::cout << "\tDist. coeff. error: " << dist_err_median << " (median), " << dist_err_mean << " (mean)" << std::endl;
    return br;
}


int main(int argc, char *argv[]) {
    /* Timing experiments */
    int nbr_iter = 1e3;
    double point_noise = 0.0;
    if (argc > 1) {
        point_noise = atof(argv[1]);
    }
    if (argc > 2) {
        nbr_iter = atoi(argv[2]);
    }
    bool print_to_file = true;
    BenchmarkResults br;
    
    HomLib::ProblemConfig config;
    config.number_points = 5;
    config.point_noise = point_noise;

    /* TEST SOLVERS HERE */
    std::cout << "================================== BENCHMARKING ==================================" << std::endl;
    std::cout << "nbr_iter = " << nbr_iter << " and point_noise = " << point_noise << std::endl;
    std::cout << "\nDOUBLE (NOT EQUAL)" << std::endl;
    config.one_sided = false;
    config.equal = false;

    std::cout << "Kukelova et al., CVPR 2015 (5 pt)" << std::endl;
    br = benchmark_solver<SolverKukelova>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "kukelova_two_sided");

    std::cout << "Kukelova et al., CVPR 2015 (6 pt)" << std::endl;
    config.number_points = 6;
    br = benchmark_solver<SolverKukelova6pt>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "kukelova_two_sided_6pt");
   
    std::cout << "Wadenback 2025 (4.5 pt)" << std::endl;  
    config.number_points = 5;
    br = benchmark_solver<SolverWadenbackDouble>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "wadenback_two_sided");

    std::cout << "\nDOUBLE (EQUAL)" << std::endl;
    config.one_sided = false;
    config.equal = true;
    
    std::cout << "Fitzgibbon, CVPR 2001 (5 pt)" << std::endl;   
    br = benchmark_solver<SolverFitzgibbon>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "fitzgibbon_equal");

    std::cout << "Kukelova et al., CVPR 2015 (5 pt - equal)" << std::endl;
    br = benchmark_solver<SolverKukelovaEqual>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "kukelova_two_sided_equal");

    std::cout << "Kukelova et al., CVPR 2015 (6 pt - equal)" << std::endl;
    config.number_points = 6;
    br = benchmark_solver<SolverKukelova6ptEqual>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "kukelova_two_sided_6pt_equal");

    std::cout << "Wadenback 2025 (4.5 pt)" << std::endl;
    config.number_points = 5;
    br = benchmark_solver<SolverWadenbackDoubleEqual>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "wadenback_two_sided_equal");
 
    std::cout << "\nSINGLE" << std::endl;
    config.one_sided = true;  // Does not matter
    config.equal = true;
    
    std::cout << "Nakano 2025, ICPR 2025 (4.5 pt)" << std::endl;
    br = benchmark_solver<SolverNakano>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "nakano_one_sided");
    
    std::cout << "Fitzgibbon, CVPR 2001 (4.5 pt)" << std::endl;   
    br = benchmark_solver<SolverFitzgibbonSingle>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "fitzgibbon_one_sided");

    std::cout << "Wadenback 2025 (4.5 pt)" << std::endl;
    br = benchmark_solver<SolverWadenbackOne>(config, nbr_iter);
    if (print_to_file)
        print_files(br, point_noise, "wadenback_one_sided");
    
    std::cout << "\n\nNOTE: Errors are log10" << std::endl;

    std::cout << "\n\nTesting nonlinear refinement - one sided" << std::endl;
    config.number_points = 12;
    HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);
    
    std::cout << "Number of points: " << inst.x1.size() << std::endl;
    
    HomLib::PoseData pd;
    pd.homography = inst.posedata.homography;
    pd.distortion_parameter = inst.posedata.distortion_parameter2;

    pd.homography(0,0) += 0.001;
    pd.homography(1,0) -= 0.001;
    pd.homography(2,0) += 0.002;
    pd.homography(0,1) += 0.001;
    pd.homography(1,1) += 0.001;
    pd.homography(2,1) -= 0.001;
    pd.homography(0,2) -= 0.0001;
    pd.homography(1,2) += 0.002;
    pd.homography(2,2) += 0.001;
    pd.distortion_parameter -= 0.0005;

    double error_before = inst.hom_error(pd.homography) + inst.dist_error(pd.distortion_parameter);
    std::cout << "before error: " << error_before << std::endl;
    HomLib::refinement_onesided(inst.x1, inst.x2, pd);
    double error_after = inst.hom_error(pd.homography) + inst.dist_error(pd.distortion_parameter);
    std::cout << "before after: " << error_after << std::endl;

    std::cout << "\n\nTesting nonlinear refinement - two sided equal" << std::endl;
    config.number_points = 12;
    config.one_sided = false;
    config.equal = true;
    inst = HomLib::generate_problem_instance(config);

    std::cout << "Number of points: " << inst.x1.size() << std::endl;

    pd.homography = inst.posedata.homography;
    pd.distortion_parameter = inst.posedata.distortion_parameter2;

    pd.homography(0,0) += 0.001;
    pd.homography(1,0) -= 0.001;
    pd.homography(2,0) += 0.002;
    pd.homography(0,1) += 0.001;
    pd.homography(1,1) += 0.001;
    pd.homography(2,1) -= 0.001;
    pd.homography(0,2) -= 0.0001;
    pd.homography(1,2) += 0.002;
    pd.homography(2,2) += 0.001;
    pd.distortion_parameter -= 0.0005;

    error_before = inst.hom_error(pd.homography) + inst.dist_error(pd.distortion_parameter);
    std::cout << "before error: " << error_before << std::endl;
    HomLib::refinement_twosided_equal(inst.x1, inst.x2, pd);
    error_after = inst.hom_error(pd.homography) + inst.dist_error(pd.distortion_parameter);
    std::cout << "before after: " << error_after << std::endl;

    std::cout << "\n\nTesting nonlinear refinement - two sided" << std::endl;
    config.number_points = 12;
    config.one_sided = false;
    config.equal = false;
    inst = HomLib::generate_problem_instance(config);

    std::cout << "Number of points: " << inst.x1.size() << std::endl;
    
    pd.homography = inst.posedata.homography;
    pd.distortion_parameter = inst.posedata.distortion_parameter;
    pd.distortion_parameter2 = inst.posedata.distortion_parameter2;

    pd.homography(0,0) += 0.001;
    pd.homography(1,0) -= 0.001;
    pd.homography(2,0) += 0.002;
    pd.homography(0,1) += 0.001;
    pd.homography(1,1) += 0.001;
    pd.homography(2,1) -= 0.001;
    pd.homography(0,2) -= 0.0001;
    pd.homography(1,2) += 0.002;
    pd.homography(2,2) += 0.001;
    pd.distortion_parameter -= 0.0005;
    pd.distortion_parameter2 += 0.0005;

    error_before = inst.hom_error(pd.homography) + inst.dist_error(pd.distortion_parameter, pd.distortion_parameter2);
    std::cout << "before error: " << error_before << std::endl;
    HomLib::refinement_twosided(inst.x1, inst.x2, pd);
    error_after = inst.hom_error(pd.homography) + inst.dist_error(pd.distortion_parameter, pd.distortion_parameter2);
    std::cout << "before after: " << error_after << std::endl;

    return 0;
}
