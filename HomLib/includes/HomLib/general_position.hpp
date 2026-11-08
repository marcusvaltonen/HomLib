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

#ifndef INCLUDES_HOMLIB_GENERAL_POSITION_HPP_
#define INCLUDES_HOMLIB_GENERAL_POSITION_HPP_

// Checks on the four points that a 4-point homography uses. Four points x[0..3] and y[0..3]
// determine a homography only if they are in general position: no three of them collinear, in
// either image. And for a plane seen by two cameras, every point is in front of both cameras, so
// the homography keeps the orientation of every triangle of the points, or flips all of them
// (when det H < 0). A mix means that some points are on the other side of the vanishing line of
// H: no physical homography maps them. Both tests use the signed areas of the triangles.

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

namespace HomLib {

// Twice the signed area of the triangle (a, b, c): positive if counterclockwise
inline double signed_area(const Eigen::Vector2d &a, const Eigen::Vector2d &b,
                          const Eigen::Vector2d &c) {
    return (b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0));
}

// A triangle counts as flat when its area is below this fraction of the largest triangle of the
// same points: this only catches exact collinearity, up to rounding
constexpr double kCollinearTolerance = 1e-10;

// Whether the four points x[0..3] and y[0..3] are a valid sample for a 4-point homography: they
// determine a homography with every point in front of both cameras. That is, no three of the points
// are collinear in either image (the points are in general position), and the four triangles all
// keep their orientation or all flip it. False for NaN points.
//
// Why the orientations: if H maps x_i to y_i, the area of each triangle (i, j, k) is multiplied by
// det(H) / (w_i w_j w_k), where w_i is the third coordinate of H x_i. So the four triangles all
// keep or all flip their orientation exactly when the four w_i have the same sign: the points are
// on the same side of the vanishing line of H, as the points of a plane seen by two cameras are.
inline bool is_valid_sample(const Eigen::Vector2d *x, const Eigen::Vector2d *y) {
    const double ax[4] = {signed_area(x[0], x[1], x[2]), signed_area(x[0], x[1], x[3]),
                          signed_area(x[0], x[2], x[3]), signed_area(x[1], x[2], x[3])};
    const double ay[4] = {signed_area(y[0], y[1], y[2]), signed_area(y[0], y[1], y[3]),
                          signed_area(y[0], y[2], y[3]), signed_area(y[1], y[2], y[3])};
    const double mx =
        std::max({std::abs(ax[0]), std::abs(ax[1]), std::abs(ax[2]), std::abs(ax[3])});
    const double my =
        std::max({std::abs(ay[0]), std::abs(ay[1]), std::abs(ay[2]), std::abs(ay[3])});
    bool ok = mx > 0 && my > 0;
    for (int t = 0; t < 4; t++) {
        // Check collinearity: the triangle isn't flat, in either image
        ok = ok && std::abs(ax[t]) > kCollinearTolerance * mx;
        ok = ok && std::abs(ay[t]) > kCollinearTolerance * my;
        // Check orientation: the four triangles all keep their orientation, or all flip it
        ok = ok && (ax[t] * ay[t] > 0) == (ax[0] * ay[0] > 0);
    }
    return ok;
}

// Chooses 4 of the 5 points x[0..4] and y[0..4] for a 4-point homography, and moves them to
// x[0..3] and y[0..3]: the 4 points whose most collinear triangle is the least collinear, in both
// images. Returns whether they are a valid sample (is_valid_sample).
inline bool select_valid_sample(Eigen::Vector2d *x, Eigen::Vector2d *y) {
    // The 10 triangles of the 5 points
    static constexpr int tri[10][3] = {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {0, 2, 3}, {0, 2, 4},
                                       {0, 3, 4}, {1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}};
    double ax[10], ay[10], mx = 0, my = 0;
    for (int t = 0; t < 10; t++) {
        ax[t] = std::abs(signed_area(x[tri[t][0]], x[tri[t][1]], x[tri[t][2]]));
        ay[t] = std::abs(signed_area(y[tri[t][0]], y[tri[t][1]], y[tri[t][2]]));
        mx = std::max(mx, ax[t]);
        my = std::max(my, ay[t]);
    }
    // All the triangles flat in an image (e.g. 5 identical points): no valid sample, and no 0 / 0
    // below
    if (!(mx > 0 && my > 0))
        return false;
    // score[d]: the most collinear triangle (relative area, the smaller in both images) of the 4
    // points other than d, i.e. of the triangles without d
    double score[5] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};
    for (int t = 0; t < 10; t++) {
        const double s = std::min(ax[t] / mx, ay[t] / my);
        for (int d = 0; d < 5; d++)
            if (d != tri[t][0] && d != tri[t][1] && d != tri[t][2])
                score[d] = std::min(score[d], s);
    }
    int drop = 4;   // the first 4 points, unless another choice is better
    for (int d = 0; d < 4; d++)
        if (score[d] > score[drop])
            drop = d;
    std::swap(x[drop], x[4]);
    std::swap(y[drop], y[4]);
    return is_valid_sample(x, y);
}

}  // namespace HomLib

#endif  // INCLUDES_HOMLIB_GENERAL_POSITION_HPP_
