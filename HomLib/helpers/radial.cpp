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
    Eigen::Vector2d radialdistort(const Eigen::Vector2d &x, double kappa) {
        // We expect inhomogenous input data
        // The inverse of radialundistort, ru = rd / (1 + kappa rd^2): rd is the root of
        // kappa ru rd^2 - rd + ru = 0 that tends to ru as kappa -> 0,
        // rd = (1 - sqrt(1 - 4 kappa ru^2)) / (2 kappa ru) = 2 ru / (1 + sqrt(1 - 4 kappa ru^2)).
        // The second form has no cancellation, and no division by kappa or ru. NaN when
        // 4 kappa ru^2 > 1: no real solution.
        return (2.0 / (1.0 + std::sqrt(1.0 - 4.0 * kappa * x.squaredNorm()))) * x;
    }
    
    bool try_radialdistort(const Eigen::Vector2d &x, double kappa, Eigen::Vector2d &out) {
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

    std::vector<Eigen::Vector2d> radialdistort(const std::vector<Eigen::Vector2d> &x, double kappa) {
        // We expect inhomogenous input data
        std::vector<Eigen::Vector2d> out;
        for (size_t i = 0; i < x.size(); i++) {
            Eigen::Vector2d y = HomLib::radialdistort(x[i], kappa);
            out.push_back(y);
        }
        return out;
    }

    Eigen::Vector2d radialundistort(const Eigen::Vector2d &x, double kappa) {
        // Compute undistorted coordinates
        Eigen::Vector2d y = x / (1 + kappa * x.squaredNorm());
        return y;
    }
    
    std::vector<Eigen::Vector2d> radialundistort(const std::vector<Eigen::Vector2d> &x, double kappa) {
        std::vector<Eigen::Vector2d> out;
        for (size_t i = 0; i < x.size(); i++) {
            Eigen::Vector2d y = HomLib::radialundistort(x[i], kappa);
            out.push_back(y);
        }
        return out;
    }
    
    void radialdistort(const std::vector<Eigen::Vector2d> &xu, std::vector<Eigen::Vector2d>* xd, double kappa) {
        std::vector<Eigen::Vector2d> tmp = radialdistort(xu, kappa);
        xd->clear();
        xd->reserve(tmp.size());
        for (size_t i = 0; i < tmp.size(); i++) {
            xd->push_back(tmp[i]);
        }
    }
    void radialundistort(const std::vector<Eigen::Vector2d> &xd, std::vector<Eigen::Vector2d>* xu, double kappa) {
        std::vector<Eigen::Vector2d> tmp = radialundistort(xd, kappa);
        xu->clear();
        xu->reserve(tmp.size());
        for (size_t i = 0; i < tmp.size(); i++) {
            xu->push_back(tmp[i]);
        }
    }
}  // namespace HomLib
