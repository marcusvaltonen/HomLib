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

#ifndef SRC_HELPERS_INVERT_HOMOGRAPHY_HPP_
#define SRC_HELPERS_INVERT_HOMOGRAPHY_HPP_

#include <Eigen/Dense>

namespace HomLib {
  // H_inv = H^-1; false (H_inv not finite) if H is singular or the inverse overflows: the solvers
  // that compute H^-1 return no solution for such an H
  inline bool invert_homography(const Eigen::Matrix3d &H, Eigen::Matrix3d *H_inv) {
    bool invertible;
    H.computeInverseWithCheck(*H_inv, invertible, 0.0);
    return invertible && H_inv->allFinite();
  }
}

#endif  // SRC_HELPERS_INVERT_HOMOGRAPHY_HPP_
