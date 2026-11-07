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

TEST_CASE("Helpers - radial distortion: the overloads with an output vector") {
    const std::vector<Eigen::Vector2d> x = {Eigen::Vector2d(0.1, -0.2), Eigen::Vector2d(0.5, 0.3),
                                            Eigen::Vector2d(-0.4, 0.6)};
    const double k = -0.2;
    std::vector<Eigen::Vector2d> xu = {Eigen::Vector2d(9, 9)}, xd;   // replaced, not appended to
    HomLib::radialundistort(x, &xu, k);
    REQUIRE(xu.size() == x.size());
    for (size_t i = 0; i < x.size(); i++)
        REQUIRE(xu[i].isApprox(HomLib::radialundistort(x[i], k)));
    HomLib::radialdistort(xu, &xd, k);
    REQUIRE(xd.size() == x.size());
    for (size_t i = 0; i < x.size(); i++)
        REQUIRE(xd[i].isApprox(x[i], 1e-12));
}
