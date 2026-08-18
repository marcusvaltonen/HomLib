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
        const std::vector<Eigen::Matrix2d> &Aff
    ) {

        double weight = 1.0;
        Eigen::Matrix<double, Eigen::Dynamic, 9> coefficients(6 * x.size(), 9);
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

            const double
	            kMinusWeightTimesX1 = -weight * x1,
	            kMinusWeightTimesY1 = -weight * y1,
	            kWeightTimesX2 = weight * x2,
	            kWeightTimesY2 = weight * y2;

            coefficients(rowIdx, 0) = kMinusWeightTimesX1;
            coefficients(rowIdx, 1) = kMinusWeightTimesY1;
            coefficients(rowIdx, 2) = -weight;
            coefficients(rowIdx, 3) = 0;
            coefficients(rowIdx, 4) = 0;
            coefficients(rowIdx, 5) = 0;
            coefficients(rowIdx, 6) = kWeightTimesX2 * x1;
            coefficients(rowIdx, 7) = kWeightTimesX2 * y1;
            coefficients(rowIdx, 8) = kWeightTimesX2;
            ++rowIdx;

            coefficients(rowIdx, 0) = 0;
            coefficients(rowIdx, 1) = 0;
            coefficients(rowIdx, 2) = 0;
            coefficients(rowIdx, 3) = kMinusWeightTimesX1;
            coefficients(rowIdx, 4) = kMinusWeightTimesY1;
            coefficients(rowIdx, 5) = -weight;
            coefficients(rowIdx, 6) = kWeightTimesY2 * x1;
            coefficients(rowIdx, 7) = kWeightTimesY2 * y1;
            coefficients(rowIdx, 8) = kWeightTimesY2;
            ++rowIdx;

            // If the minimal case is considered, we 
            // do not need all constraints to estimate 
            // the homography.
            //if (i == 1) {
	        //    break;
            //}
            
            // NOTE(MARCUS): Degenerates without.... don't know why

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
        
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(coefficients,Eigen::ComputeThinV);

        const Eigen::Matrix<double, 9, 1> &h = svd.matrixV().rightCols<1>();

        Eigen::Matrix3d H;
        H << h(0), h(1), h(2),
             h(3), h(4), h(5),
             h(6), h(7), h(8);
        H /= H(2,2);

        std::vector<HomLib::PoseData> output;
        HomLib::PoseData pd;
        pd.homography = H;
        pd.distortion_parameter = 0.0;
        pd.distortion_parameter2 = 0.0;
        output.push_back(pd);

        return output;
    }
}  // namespace BarathVISAPP2016
}  // namespace HomLib
