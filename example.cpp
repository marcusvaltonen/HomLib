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
#include <chrono>  // NOLINT [build/c++11]
#include <iostream>
#include "get_fitzgibbon_cvpr_2001.hpp"
#include "get_kukelova_cvpr_2015.hpp"
#include "get_valtonenornhag_arxiv_2020a.hpp"
#include "get_valtonenornhag_arxiv_2020b.hpp"
#include "problem_instance.hpp"
#include "generate_problem_instance.hpp"
#include "posedata.hpp"

int main() {
    /* Timing experiments */
    int nbr_iter = 1e2;
    Eigen::MatrixXcd out;
    int N = 5;
    HomLib::ProblemInstance inst;

    /* TEST SOLVERS HERE */
    std::cout << "--- COMPLETE SOLVERS, INCLUDING PRE- and POST-PROCESSING ---" << std::endl;
    HomLib::PoseData posedata;
    N = 2;

    // fHf complete
    std::vector<HomLib::PoseData> posedata_fHf;
    auto start = std::chrono::steady_clock::now();
    for (int i=0; i < nbr_iter; i++) {
        inst = HomLib::generate_problem_instance(N);
        posedata_fHf = HomLib::ValtonenOrnhagArxiv2020B::get_fHf(inst.x1, inst.x2, inst.R1g, inst.R2g);
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "(fHf complete) Elapsed time in microseconds : "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << std::endl;

    // frHfr elim complete
    N = 3;
    start = std::chrono::steady_clock::now();
    for (int i=0; i < nbr_iter; i++) {
        inst = HomLib::generate_problem_instance(N);
        posedata = HomLib::ValtonenOrnhagArxiv2020B::get_frHfr(inst.x1, inst.x2, inst.R1g, inst.R2g);
    }
    end = std::chrono::steady_clock::now();
    std::cout << "(frHfr elim complete) Elapsed time in microseconds : "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << std::endl;

    // floor_fHf complete
    start = std::chrono::steady_clock::now();
    for (int i=0; i < nbr_iter; i++) {
        inst = HomLib::generate_problem_instance(N);
        posedata = HomLib::ValtonenOrnhagArxiv2020A::get_fHf(inst.x1, inst.x2, inst.R1g, inst.R2g);
    }
    end = std::chrono::steady_clock::now();
    std::cout << "(floor fHf complete) Elapsed time in microseconds : "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << std::endl;

    // Kukelova Hlam1lam2 complete
    std::vector<HomLib::PoseData> posedata2;
    start = std::chrono::steady_clock::now();
    N = 5;
    for (int i=0; i < nbr_iter; i++) {
        inst = HomLib::generate_problem_instance(N);
        posedata2 = HomLib::KukelovaCVPR2015::get(inst.x1, inst.x2);
    }
    end = std::chrono::steady_clock::now();
    std::cout << "(Kukelova Hlam1lam2 complete) Elapsed time in microseconds : "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << std::endl;

    // Fitzgibbon
    start = std::chrono::steady_clock::now();
    for (int i=0; i < nbr_iter; i++) {
        inst = HomLib::generate_problem_instance(N);
        posedata = HomLib::FitzgibbonCVPR2001::get(inst.x1, inst.x2);
    }
    end = std::chrono::steady_clock::now();
    std::cout << "(Fitzgibbon complete) Elapsed time in microseconds : "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / nbr_iter
        << " µs" << std::endl;

    return 0;
}
