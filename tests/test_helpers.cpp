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
#include <catch2/catch_test_macros.hpp>
#include <vector>
#include <cmath>
#include <limits>
#include "radial.hpp"
#include "roots.hpp"

TEST_CASE("Helpers - roots") {
    Eigen::VectorXd x(5);
    x << 1, 2, 3, 4, 5;
    Eigen::MatrixXd p2(2, 5);

    Eigen::VectorXcd w = HomLib::roots(x);

    typedef std::complex<double> cd;
    Eigen::VectorXcd expected(4);
    expected << cd(0.287815479557648, 1.416093080171908),
                cd(0.287815479557648, -1.416093080171908),
                cd(-1.287815479557648, 0.857896758328491),
                cd(-1.287815479557648, -0.857896758328491);

    double tol = 1e-12;
    REQUIRE(w.isApprox(expected, tol));
}

TEST_CASE("Helpers - radial distortion: try_radialdistort inverts try_radialundistort, also at the centre") {
    for (double k : {-0.5, -0.1, 0.0, 0.1, 0.2}) {
        Eigen::Vector2d out;
        REQUIRE(HomLib::try_radialdistort(Eigen::Vector2d(0, 0), k, out));
        REQUIRE(out == Eigen::Vector2d(0, 0));
        for (double r : {1e-9, 1e-4, 0.3, 1.0}) {
            const Eigen::Vector2d xd(0.6 * r, -0.8 * r);
            Eigen::Vector2d xu, back;
            REQUIRE(HomLib::try_radialundistort(xd, k, xu));
            REQUIRE(HomLib::try_radialdistort(xu, k, back));
            REQUIRE((back - xd).norm() <= 1e-15 * r);
        }
    }
}

TEST_CASE("Helpers - radial distortion: the versions that report failure") {
    const Eigen::Vector2d unchanged(7, 7);
    for (double k : {-0.5, 0.0, 0.2}) {
        for (const Eigen::Vector2d x : {Eigen::Vector2d(0, 0), Eigen::Vector2d(0.3, -0.4), Eigen::Vector2d(1.0, 0.0)}) {
            Eigen::Vector2d out = unchanged;
            REQUIRE(HomLib::try_radialundistort(x, k, out));
            REQUIRE(out == x / (1 + k * x.squaredNorm()));
            out = unchanged;
            REQUIRE(HomLib::try_radialdistort(x, k, out));
            REQUIRE((out / (1 + k * out.squaredNorm()) - x).norm() <= 1e-15);
        }
    }
    Eigen::Vector2d out = unchanged;
    // At the radius 1 / sqrt(-k) where the model folds: 1 + k |x|^2 = 0. Beyond it, the
    // undistorted point is finite (the same projective point as (-x, -w)).
    REQUIRE_FALSE(HomLib::try_radialundistort(Eigen::Vector2d(1.0, 0.0), -1.0, out));
    REQUIRE_FALSE(HomLib::try_radialundistort(Eigen::Vector2d(4.0, 0.0), -0.0625, out));
    Eigen::Vector2d beyond;
    REQUIRE(HomLib::try_radialundistort(Eigen::Vector2d(2.0, 0.0), -0.5, beyond));
    REQUIRE(beyond == Eigen::Vector2d(-2.0, 0.0));
    // No distorted point maps to x: 4 k |x|^2 > 1
    REQUIRE_FALSE(HomLib::try_radialdistort(Eigen::Vector2d(2.0, 0.0), 0.1, out));
    // |x|^2 overflows, or a NaN
    REQUIRE_FALSE(HomLib::try_radialundistort(Eigen::Vector2d(1e200, 0.0), 0.0, out));
    REQUIRE_FALSE(HomLib::try_radialdistort(Eigen::Vector2d(1e200, 0.0), 0.0, out));
    REQUIRE_FALSE(HomLib::try_radialdistort(Eigen::Vector2d(0.5, 0.0), std::numeric_limits<double>::quiet_NaN(), out));
    REQUIRE(out == unchanged);
}
