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


#include <math.h>
#include <random>

#include "problem_instance.hpp"
#include "generate_problem_instance.hpp"

namespace HomLib {
    inline Eigen::Matrix3d rotx(double angle);
    inline Eigen::Matrix3d rotz(double angle);
    HomLib::ProblemInstance generate_problem_instance(int N) {
        HomLib::ProblemInstance p;

        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<double> d{0, 1};

        double tmp1, tmp2;
        tmp1 = d(gen);
        tmp2 = d(gen);

        Eigen::Matrix3d R1g, R2g;
        Eigen::MatrixXd x1(3, N), x2(3, N);
        R1g = rotx(tmp1) * rotz(tmp2);

        tmp1 = d(gen);
        tmp2 = d(gen);

        R2g = rotx(tmp1) * rotz(tmp2);

        // FIXME: For speed tests this does not matter.. but it is not on a plane
        x1 = Eigen::MatrixXd::Random(2, N);
        x2 = Eigen::MatrixXd::Random(2, N);

        p.R1g = R1g;
        p.R2g = R2g;
        p.x1 = x1;
        p.x2 = x2;

        return p;
    }

    inline Eigen::Matrix3d rotx(double angle) {
        Eigen::Matrix3d R;
        R << 1,           0,            0,
             0, std::cos(angle), -std::sin(angle),
             0, std::sin(angle),  std::cos(angle);
        return R;
    }

    inline Eigen::Matrix3d rotz(double angle) {
        Eigen::Matrix3d R;
        R << std::cos(angle), -std::sin(angle), 0,
             std::sin(angle),  std::cos(angle), 0,
             0          ,            0, 1;
        return R;
    }
}  // namespace HomLib
