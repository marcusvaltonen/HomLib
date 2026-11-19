// Copyright (c) 2020, Viktor Larsson
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of the copyright holder nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copied from PoseLib (https://github.com/PoseLib/PoseLib), PoseLib/misc/univariate.cc at commit
// a69263d5d824b5de11e48756f33f6e8ebb231f0c: solve_cubic_real only, unchanged except the namespace
// (poselib::univariate there).

#include "univariate.hpp"

#include <cmath>

namespace HomLib {
namespace univariate {

int solve_cubic_real(double c2, double c1, double c0, double roots[3]) {
    double a = c1 - c2 * c2 / 3.0;
    double b = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1) / 27.0 + c0;
    double c = b * b / 4.0 + a * a * a / 27.0;
    int n_roots;
    if (a == 0.0 && b == 0.0) {
        // Triple root at -c2/3
        roots[0] = roots[1] = roots[2] = -c2 / 3.0;
        n_roots = 3;
    } else if (c > 0) {
        c = std::sqrt(c);
        b *= -0.5;
        roots[0] = std::cbrt(b + c) + std::cbrt(b - c) - c2 / 3.0;
        n_roots = 1;
    } else {
        c = 3.0 * b / (2.0 * a) * std::sqrt(-3.0 / a);
        double d = 2.0 * std::sqrt(-a / 3.0);
        roots[0] = d * std::cos(std::acos(c) / 3.0) - c2 / 3.0;
        roots[1] = d * std::cos(std::acos(c) / 3.0 - 2.09439510239319526263557236234192) - c2 / 3.0; // 2*pi/3
        roots[2] = d * std::cos(std::acos(c) / 3.0 - 4.18879020478639052527114472468384) - c2 / 3.0; // 4*pi/3
        n_roots = 3;
    }

    // single newton iteration
    for (int i = 0; i < n_roots; ++i) {
        double x = roots[i];
        double x2 = x * x;
        double x3 = x * x2;
        double dx = -(x3 + c2 * x2 + c1 * x + c0) / (3 * x2 + 2 * c2 * x + c1);
        roots[i] += dx;
    }
    return n_roots;
}

} // namespace univariate
} // namespace HomLib
