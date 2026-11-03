#pragma once
#if defined(HOMLIB_HAS_GCRANSAC) || defined(HOMLIB_HAS_MAGSAC)

#include <vector>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

namespace HomLib {
namespace detail {

inline cv::Mat pack_point_cvmat(
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y)
{
    const int n = static_cast<int>(x.size());
    cv::Mat data(n, 4, CV_64F);
    for (int i = 0; i < n; ++i) {
        double* row = data.ptr<double>(i);
        row[0] = x[i](0); row[1] = x[i](1);
        row[2] = y[i](0); row[3] = y[i](1);
    }
    return data;
}

inline cv::Mat pack_affine_cvmat(
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Matrix2d>& A)
{
    const int n = static_cast<int>(x.size());
    cv::Mat data(n, 8, CV_64F);
    for (int i = 0; i < n; ++i) {
        double* row = data.ptr<double>(i);
        row[0] = x[i](0); row[1] = x[i](1);
        row[2] = y[i](0); row[3] = y[i](1);
        row[4] = A[i](0, 0); row[5] = A[i](0, 1);
        row[6] = A[i](1, 0); row[7] = A[i](1, 1);
    }
    return data;
}

inline cv::Mat pack_orientation_cvmat(
    const std::vector<Eigen::Vector2d>& x,
    const std::vector<Eigen::Vector2d>& y,
    const std::vector<Eigen::Vector2d>& ori)
{
    const int n = static_cast<int>(x.size());
    cv::Mat data(n, 6, CV_64F);
    for (int i = 0; i < n; ++i) {
        double* row = data.ptr<double>(i);
        row[0] = x[i](0); row[1] = x[i](1);
        row[2] = y[i](0); row[3] = y[i](1);
        row[4] = ori[i](0); row[5] = ori[i](1);
    }
    return data;
}

} // namespace detail
} // namespace HomLib

#endif // HOMLIB_HAS_GCRANSAC || HOMLIB_HAS_MAGSAC
