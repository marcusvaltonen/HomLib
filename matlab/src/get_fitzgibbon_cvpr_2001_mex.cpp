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
#include <vector>
#include "get_fitzgibbon_cvpr_2001.hpp"

#ifdef MATLAB_MEX_FILE
#include "mex.h"  // NOLINT [build/include_subdir]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("get_fitzgibbon_cvpr_2001:nrhs", "Two inputs required.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("get_fitzgibbon_cvpr_2001:nlhs", "Two outputs required.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("get_fitzgibbon_cvpr_2001:notDouble", "Input data must be type double.");
    }
    if (mxGetNumberOfElements(prhs[0]) != 10) {
        mexErrMsgIdAndTxt("get_fitzgibbon_cvpr_2001:incorrectSize1", "Incorrect input size of first argument");
    }
    if (mxGetNumberOfElements(prhs[1]) != 10) {
        mexErrMsgIdAndTxt("get_fitzgibbon_cvpr_2001:incorrectSize2", "Incorrect input size of second argument");
    }

    // Convert to expected input
    // TODO(marcusvaltonen): Cast directly to MatrixXd?
    Eigen::VectorXd x1_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[0]), 10);
    Eigen::VectorXd x2_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[1]), 10);
    Eigen::MatrixXd x1 = Eigen::Map<Eigen::MatrixXd>(x1_tmp.data(), 2, 5);
    Eigen::MatrixXd x2 = Eigen::Map<Eigen::MatrixXd>(x2_tmp.data(), 2, 5);

    // Compute output
    HomLib::PoseData posedata = HomLib::FitzgibbonCVPR2001::get(x1, x2);

    // Wrap it up to Matlab compatible output
    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* zr = mxGetPr(plhs[0]);
    for (Eigen::Index i = 0; i < posedata.homography.size(); i++) {
        zr[i] = posedata.homography(i);
    }

    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    zr = mxGetPr(plhs[1]);
    zr[0] = posedata.distortion_parameter;
}
#endif  // MATLAB_MEX_FILE
