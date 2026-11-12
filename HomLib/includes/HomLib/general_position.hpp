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

// A triangle counts as flat when its area is below this fraction of the size of the points (see
// is_valid_sample): this only catches collinearity, up to rounding
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
//
// Flat triangles: a triangle is flat when its area is below kCollinearTolerance times the size of
// the four points, the sum of the 6 squared distances between them. The size has the units of an
// area and doesn't depend on the order of the points. Unlike the largest of the four triangles, it
// doesn't vanish when all four points are collinear, so such samples are rejected. (It is also
// stricter when one point is very far from the other three: they then count as collinear once
// their distances are below about 2e-5 times that distance.)
//
// Optimizations (the same decisions, 2-2.5x faster):
// - Shared differences: everything comes from the 3 differences d_k = p_k - p_0 of each image. The
//   areas of the 3 triangles with p_0 are the cross products of the d_k, and the 4th triangle's is
//   the cross product of d_2 - d_1 and d_3 - d_1.
// - The size without the 6 distances: sum_{i<j} |p_i - p_j|^2 = 4 sum_k |d_k|^2 - |sum_k d_k|^2.
//   (Both sides are 4 times the spread of the points about their centroid, so p_0 isn't special.)
// - No branches: the checks are accumulated with &= instead of chained with &&, so the compiler
//   generates straight-line code (no jumps with GCC or clang at -O3; clang also processes both
//   images at once, with 2-lane SIMD instructions).
inline bool is_valid_sample(const Eigen::Vector2d *x, const Eigen::Vector2d *y) {
    // For the points p[0..3]: the signed areas (twice) a[0..3] of the triangles (0, 1, 2),
    // (0, 1, 3), (0, 2, 3), (1, 2, 3), and the tolerance below which a triangle is flat
    auto areas = [](const Eigen::Vector2d *p, double *a, double &tolerance) {
        const double u1 = p[1](0) - p[0](0), v1 = p[1](1) - p[0](1);
        const double u2 = p[2](0) - p[0](0), v2 = p[2](1) - p[0](1);
        const double u3 = p[3](0) - p[0](0), v3 = p[3](1) - p[0](1);
        a[0] = u1 * v2 - v1 * u2;
        a[1] = u1 * v3 - v1 * u3;
        a[2] = u2 * v3 - v2 * u3;
        a[3] = (u2 - u1) * (v3 - v1) - (v2 - v1) * (u3 - u1);
        const double su = u1 + u2 + u3, sv = v1 + v2 + v3;
        const double size = 4 * (u1 * u1 + v1 * v1 + u2 * u2 + v2 * v2 + u3 * u3 + v3 * v3) -
                            su * su - sv * sv;
        tolerance = kCollinearTolerance * size;
    };
    double ax[4], ay[4], tx, ty;
    areas(x, ax, tx);
    areas(y, ay, ty);
    // False if all the points of an image are equal, or for a NaN
    bool ok = (tx > 0) & (ty > 0);
    const bool orientation = ax[0] * ay[0] > 0;
    for (int t = 0; t < 4; t++) {
        // Check collinearity: the triangle isn't flat, in either image
        ok &= std::abs(ax[t]) > tx;
        ok &= std::abs(ay[t]) > ty;
        // Check orientation: the four triangles all keep their orientation, or all flip it
        ok &= (ax[t] * ay[t] > 0) == orientation;
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
