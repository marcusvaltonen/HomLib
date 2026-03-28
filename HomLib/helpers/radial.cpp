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
#include <math.h>  // copysign
#include <cmath>

#include "radial.hpp"

namespace HomLib {
    Eigen::Vector2d radialdistort(const Eigen::Vector2d &x, double kappa) {
        // We expect inhomogenous input data

        double ru2 = x.squaredNorm();
        double ru = std::sqrt(ru2);

        double rd;
        // Compute distorted radius
        if (kappa == 0) {
            rd = ru;
        } else {
            rd = 0.5 / kappa / ru - copysign(1.0, kappa) * std::sqrt(0.25 / std::pow(kappa, 2) / ru2 - 1.0 / kappa);
        }

        // Compute distorted coordinates
        Eigen::Vector2d y = (rd / ru) * x;

        return y;
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
        for (size_t i = 0; tmp.size(); i++) {
            xu->push_back(tmp[i]);
        }
    }
}  // namespace HomLib
