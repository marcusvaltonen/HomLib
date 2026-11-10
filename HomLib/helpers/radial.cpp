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
#include <cmath>

#include "radial.hpp"

namespace HomLib {
    bool try_radialdistort(const Eigen::Vector2d &x, double kappa, Eigen::Vector2d &out) {
        // The inverse of try_radialundistort, ru = rd / (1 + kappa rd^2): rd is the root of
        // kappa ru rd^2 - rd + ru = 0 that tends to ru as kappa -> 0,
        // rd = (1 - sqrt(1 - 4 kappa ru^2)) / (2 kappa ru) = 2 ru / (1 + sqrt(1 - 4 kappa ru^2)).
        // The second form has no cancellation, and no division by kappa or ru; no real solution
        // when 4 kappa ru^2 > 1.
        const double r2 = x.squaredNorm();
        if (!std::isfinite(r2))
            return false;
        const double d = 1.0 - 4.0 * kappa * r2;
        if (!std::isgreaterequal(d, 0.0))   // quiet comparison: no exception if d is NaN
            return false;
        out = (2.0 / (1.0 + std::sqrt(d))) * x;
        return true;
    }

    bool try_radialundistort(const Eigen::Vector2d &x, double kappa, Eigen::Vector2d &out) {
        const double r2 = x.squaredNorm();
        if (!std::isfinite(r2))
            return false;
        const double w = 1.0 + kappa * r2;
        if (!std::isfinite(w) || w == 0.0)
            return false;
        out = x / w;
        return true;
    }
}  // namespace HomLib
