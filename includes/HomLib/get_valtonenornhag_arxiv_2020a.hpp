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

#ifndef INCLUDES_HOMLIB_GET_VALTONENORNHAG_ARXIV_2020A_HPP_
#define INCLUDES_HOMLIB_GET_VALTONENORNHAG_ARXIV_2020A_HPP_

#include <Eigen/Dense>
#include "posedata.hpp"

namespace HomLib::ValtonenOrnhagArxiv2020A {
    HomLib::PoseData get_fHf(
        const Eigen::MatrixXd &x1,
        const Eigen::MatrixXd &x2,
        const Eigen::Matrix3d &R1,
        const Eigen::Matrix3d &q2);
}

#endif  // INCLUDES_HOMLIB_GET_VALTONENORNHAG_ARXIV_2020A_HPP_
