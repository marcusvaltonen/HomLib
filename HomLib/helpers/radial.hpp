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

#ifndef SRC_HELPERS_RADIAL_HPP_
#define SRC_HELPERS_RADIAL_HPP_

#include <Eigen/Dense>

namespace HomLib {
  Eigen::Vector2d radialdistort(const Eigen::Vector2d &x, double kappa);
  Eigen::Vector2d radialundistort(const Eigen::Vector2d &x, double kappa);
  std::vector<Eigen::Vector2d> radialdistort(const std::vector<Eigen::Vector2d> &x, double kappa);
  std::vector<Eigen::Vector2d> radialundistort(const std::vector<Eigen::Vector2d> &x, double kappa);
  void radialdistort(const std::vector<Eigen::Vector2d> &xu, std::vector<Eigen::Vector2d>* xd, double kappa);
  void radialundistort(const std::vector<Eigen::Vector2d> &xd, std::vector<Eigen::Vector2d>* xu, double kappa);

  // Versions that return false, without an invalid operation or a division by zero (they don't
  // trap with floating-point exceptions enabled), and leave out unchanged: try_radialundistort
  // where 1 + kappa |x|^2 = 0, try_radialdistort where 4 kappa |x|^2 > 1 (no distorted point maps
  // to x), both if |x|^2 overflows.
  bool try_radialundistort(const Eigen::Vector2d &x, double kappa, Eigen::Vector2d &out);
  bool try_radialdistort(const Eigen::Vector2d &x, double kappa, Eigen::Vector2d &out);
}

#endif  // SRC_HELPERS_RADIAL_HPP_
