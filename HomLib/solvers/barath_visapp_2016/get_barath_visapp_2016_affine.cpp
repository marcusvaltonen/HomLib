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


#include <vector>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include "posedata.hpp"

#include "radial.hpp"
#include "get_barath_visapp_2016.hpp"

namespace HomLib {
namespace BarathVISAPP2016 {
    std::vector<HomLib::PoseData> get_affine(
        const std::vector<Eigen::Vector2d> &x,
        const std::vector<Eigen::Vector2d> &y,
        const std::vector<Eigen::Matrix2d> &Aff,
        bool use_qr
    ) {

        const bool minimal = (x.size() == 2);
        const int nrows = minimal ? 8 : static_cast<int>(6 * x.size());
        Eigen::Matrix<double, Eigen::Dynamic, 9> coefficients(nrows, 9);
        size_t rowIdx = 0;
        for (size_t i = 0; i < x.size(); i++)
        {
            const double
	            x1 = x[i](0),
	            y1 = x[i](1),
	            x2 = y[i](0),
	            y2 = y[i](1),
	            a11 = Aff[i](0,0),
	            a12 = Aff[i](0,1),
	            a21 = Aff[i](1,0),
	            a22 = Aff[i](1,1);

            coefficients(rowIdx, 0) = -x1;
            coefficients(rowIdx, 1) = -y1;
            coefficients(rowIdx, 2) = -1;
            coefficients(rowIdx, 3) = 0;
            coefficients(rowIdx, 4) = 0;
            coefficients(rowIdx, 5) = 0;
            coefficients(rowIdx, 6) = x2 * x1;
            coefficients(rowIdx, 7) = x2 * y1;
            coefficients(rowIdx, 8) = x2;
            ++rowIdx;

            coefficients(rowIdx, 0) = 0;
            coefficients(rowIdx, 1) = 0;
            coefficients(rowIdx, 2) = 0;
            coefficients(rowIdx, 3) = -x1;
            coefficients(rowIdx, 4) = -y1;
            coefficients(rowIdx, 5) = -1;
            coefficients(rowIdx, 6) = y2 * x1;
            coefficients(rowIdx, 7) = y2 * y1;
            coefficients(rowIdx, 8) = y2;
            ++rowIdx;

            // a11 row
            coefficients(rowIdx, 0) = -1;
            coefficients(rowIdx, 1) = 0;
            coefficients(rowIdx, 2) = 0;
            coefficients(rowIdx, 3) = 0;
            coefficients(rowIdx, 4) = 0;
            coefficients(rowIdx, 5) = 0;
            coefficients(rowIdx, 6) = x2 + a11 * x1;
            coefficients(rowIdx, 7) = a11 * y1;
            coefficients(rowIdx, 8) = a11;
            ++rowIdx;

            // Minimal: AC2 contributes only 1 affine row (a11)
            if (minimal && i >= 1)
                continue;

            // a12 row
            coefficients(rowIdx, 0) = 0;
            coefficients(rowIdx, 1) = -1;
            coefficients(rowIdx, 2) = 0;
            coefficients(rowIdx, 3) = 0;
            coefficients(rowIdx, 4) = 0;
            coefficients(rowIdx, 5) = 0;
            coefficients(rowIdx, 6) = a12 * x1;
            coefficients(rowIdx, 7) = x2 + a12 * y1;
            coefficients(rowIdx, 8) = a12;
            ++rowIdx;

            // a21 row
            coefficients(rowIdx, 0) = 0;
            coefficients(rowIdx, 1) = 0;
            coefficients(rowIdx, 2) = 0;
            coefficients(rowIdx, 3) = -1;
            coefficients(rowIdx, 4) = 0;
            coefficients(rowIdx, 5) = 0;
            coefficients(rowIdx, 6) = y2 + a21 * x1;
            coefficients(rowIdx, 7) = a21 * y1;
            coefficients(rowIdx, 8) = a21;
            ++rowIdx;

            // a22 row — skipped for AC1 in minimal case
            if (minimal && i == 0)
                continue;

            coefficients(rowIdx, 0) = 0;
            coefficients(rowIdx, 1) = 0;
            coefficients(rowIdx, 2) = 0;
            coefficients(rowIdx, 3) = 0;
            coefficients(rowIdx, 4) = -1;
            coefficients(rowIdx, 5) = 0;
            coefficients(rowIdx, 6) = a22 * x1;
            coefficients(rowIdx, 7) = y2 + a22 * y1;
            coefficients(rowIdx, 8) = a22;
            ++rowIdx;
        }
        
        Eigen::Matrix<double, 9, 1> h;
        if (use_qr && x.size() == 2) {
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(coefficients.transpose());
            h = qr.householderQ() * Eigen::Matrix<double, 9, 1>::Unit(8);
        } else {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(coefficients, Eigen::ComputeFullV);
            h = svd.matrixV().col(8);
        }

        Eigen::Matrix3d H;
        H << h(0), h(1), h(2),
             h(3), h(4), h(5),
             h(6), h(7), h(8);
        std::vector<HomLib::PoseData> output;
        if (H(2,2) != 0.0)   // otherwise the unit-norm h (H maps the origin to infinity)
            H /= H(2,2);
        if (!H.allFinite())
            return output;

        HomLib::PoseData pd;
        pd.homography = H;
        pd.distortion_parameter = 0.0;
        pd.distortion_parameter2 = 0.0;
        output.push_back(pd);

        return output;
    }
}  // namespace BarathVISAPP2016
}  // namespace HomLib
