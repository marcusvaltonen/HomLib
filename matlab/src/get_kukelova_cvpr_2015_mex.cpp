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
#include "get_kukelova_cvpr_2015.hpp"

#ifdef MATLAB_MEX_FILE
#include "mex.h"  // NOLINT [build/include_subdir]
#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("get_kukelova_cvpr_2015:nrhs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("get_kukelova_cvpr_2015:nlhs", "One output required.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("get_kukelova_cvpr_2015:notDouble", "Input data must be type double.");
    }
    if (mxGetNumberOfElements(prhs[0]) != 10) {
        mexErrMsgIdAndTxt("get_kukelova_cvpr_2015:incorrectSize1", "Incorrect input size of first argument");
    }
    if (mxGetNumberOfElements(prhs[1]) != 10) {
        mexErrMsgIdAndTxt("get_kukelova_cvpr_2015:incorrectSize2", "Incorrect input size of second argument");
    }

    // Convert to expected input
    // TODO(marcusvaltonen): Cast directly to MatrixXd?
    Eigen::VectorXd x1_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[0]), 10);
    Eigen::VectorXd x2_tmp = Eigen::Map<Eigen::VectorXd>(mxGetPr(prhs[1]), 10);
    Eigen::MatrixXd x1 = Eigen::Map<Eigen::MatrixXd>(x1_tmp.data(), 2, 5);
    Eigen::MatrixXd x2 = Eigen::Map<Eigen::MatrixXd>(x2_tmp.data(), 2, 5);

    // Compute output
    std::vector<HomLib::PoseData> posedata = HomLib::KukelovaCVPR2015::get(x1, x2);

    // Wrap it up to Matlab compatible output
    std::size_t NUMBER_OF_STRUCTS = posedata.size();
    const char *field_names[] = {"H", "lam1", "lam2"};
    mwSize dims[2] = {1, NUMBER_OF_STRUCTS };
    int H_field, lam1_field, lam2_field;
    mwIndex i;

    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);

    H_field = mxGetFieldNumber(plhs[0], "H");
    lam1_field = mxGetFieldNumber(plhs[0], "lam1");
    lam2_field = mxGetFieldNumber(plhs[0], "lam2");

    double* zr;
    for (i = 0; i < NUMBER_OF_STRUCTS; i++) {
        mxArray *field_value;

        // Create H
        field_value = mxCreateDoubleMatrix(3, 3, mxREAL);
        zr = mxGetPr(field_value);
        for (Eigen::Index j = 0; j < posedata[i].homography.size(); j++) {
            zr[j] = posedata[i].homography(j);
        }
        mxSetFieldByNumber(plhs[0], i, H_field, field_value);

        // Create lam1
        field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
        zr = mxGetPr(field_value);
        zr[0] = posedata[i].distortion_parameter;
        mxSetFieldByNumber(plhs[0], i, lam1_field, field_value);

        // Create lam2
        field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
        zr = mxGetPr(field_value);
        zr[0] = posedata[i].distortion_parameter2;
        mxSetFieldByNumber(plhs[0], i, lam2_field, field_value);
    }
}
#endif  // MATLAB_MEX_FILE
