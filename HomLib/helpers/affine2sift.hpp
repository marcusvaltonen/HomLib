#ifndef SRC_HELPERS_AFFINE2SIFT_HPP_
#define SRC_HELPERS_AFFINE2SIFT_HPP_

#include <vector>
#include <numeric>
#include <Eigen/Dense>

using namespace Eigen;

// Converts an affine correspondence to a scale-and-orientation correspondence.
// There are eight possible solutions, but this arbitrarily returns the first one.
void affine2sift(const Eigen::Matrix2d &A, // input affine transformation matrix
                 double &s_ref, double &c_ref, // sine and cosine of feature orientation in reference image
                 double &s_query, double &c_query, // sine and cosine of feature orientation in query image
                 double &q // ratio of feature scales (scale in query image / scale in reference image)
);

#endif  // SRC_HELPERS_AFFINE2SIFT_HPP_