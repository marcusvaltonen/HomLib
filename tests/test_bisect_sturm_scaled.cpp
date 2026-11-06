// Copyright (c) 2026 Frédéric Devernay
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

// Tests of bisect_sturm_scaled (wadenback_common.hpp), the root finder of the two-sided
// Wadenbäck solvers. Each test builds a polynomial from known roots.

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "wadenback_common.hpp"

namespace {

// Returns the ascending coefficients of lead * prod (x - r) over the real roots r, times the
// quadratic factors {a, b, c} = a + b x + c x^2. The quadratic factors must have no real roots.
std::vector<double> poly_from_roots(const std::vector<double> &real,
                                    const std::vector<std::vector<double>> &quadratics,
                                    double lead = 1.0) {
    std::vector<double> p = {lead};
    auto multiply = [&p](const std::vector<double> &f) {
        std::vector<double> q(p.size() + f.size() - 1, 0.0);
        for (size_t i = 0; i < p.size(); i++)
            for (size_t j = 0; j < f.size(); j++)
                q[i + j] += p[i] * f[j];
        p = q;
    };
    for (double r : real)
        multiply({-r, 1.0});
    for (const std::vector<double> &f : quadratics)
        multiply(f);
    return p;
}

// Returns the roots that bisect_sturm_scaled<N> finds, sorted
template <int N>
std::vector<double> scaled_roots(const std::vector<double> &p) {
    REQUIRE(p.size() == N + 1);
    double roots[N];
    const int n = HomLib::Wadenback3DV2026::bisect_sturm_scaled<N>(p.data(), roots);
    std::vector<double> found(roots, roots + n);
    std::sort(found.begin(), found.end());
    return found;
}

// Checks that found holds the expected roots and nothing else, each within
// rel_tol * max(|root|, 1e-3)
void check_roots(std::vector<double> expected, const std::vector<double> &found, double rel_tol) {
    std::sort(expected.begin(), expected.end());
    REQUIRE(found.size() == expected.size());
    for (size_t i = 0; i < found.size(); i++) {
        CAPTURE(i, expected[i], found[i]);
        const double tol = rel_tol * std::max(std::abs(expected[i]), 1e-3);
        REQUIRE(std::abs(found[i] - expected[i]) <= tol);
    }
}

}  // namespace

TEST_CASE("bisect_sturm_scaled - finds small real roots, 0 included") {
    // Degree 9: these small real roots times x^2 - 2x + 2 (roots 1 +/- i). PoseLib's bisect_sturm
    // alone finds them with errors up to about 4e-4. A problem without distortion in
    // test_wadenback_3dv_2026.cpp tests the root at 0 that bisect_sturm loses.
    const std::vector<double> real = {-0.05, -0.003, 0.0, 1e-4, 0.002, 0.01, 0.04};
    check_roots(real, scaled_roots<9>(poly_from_roots(real, {{2.0, -2.0, 1.0}})), 1e-9);
}

TEST_CASE("bisect_sturm_scaled - finds each root of a cluster of close roots") {
    // Three roots within 2e-3 of each other, as in some polynomials of noisy minimal samples. Near
    // these roots, |p| stays small over a wide interval. With its default tolerance, bisect_sturm
    // stops early there (errors up to about 4e-4). Rounding the coefficients alone moves these
    // roots by up to about 2e-10.
    const std::vector<double> real = {-0.544, -0.4442, -0.443, -0.4422, -0.3236, 2.8};
    check_roots(real, scaled_roots<6>(poly_from_roots(real, {}, 0.34)), 1e-7);
}

TEST_CASE("bisect_sturm_scaled - returns no root when there is none") {
    const std::vector<double> p =
        poly_from_roots({}, {{1.0, 0.0, 1.0}, {0.01, 0.1, 1.0}, {4.0, -1.0, 1.0}});
    REQUIRE(scaled_roots<6>(p).empty());
}

TEST_CASE("bisect_sturm_scaled - finds the roots of random polynomials") {
    // 2000 polynomials of degree 6 and 9. The real roots have random signs and magnitudes from 1e-3
    // to 1e3 (log-uniform), and are at least 5% apart. Random quadratic factors without real roots
    // fill the remaining degree. The uniform values come from the raw output of a fixed generator,
    // so every platform builds the same polynomials (up to std::pow's last bit). The test checks
    // the accuracy over a wide range of scales: without the scaling, many roots are off by more
    // than 1e-8.
    std::mt19937 rng(20261105);
    auto uniform = [&rng]() { return rng() / 4294967296.0; };   // in [0, 1)
    for (int t = 0; t < 2000; t++) {
        const int degree = t % 2 ? 9 : 6;
        const int n_real = degree - 2 * static_cast<int>(uniform() * (degree / 2 + 1));
        std::vector<double> real;
        while (static_cast<int>(real.size()) < n_real) {
            const double sign = uniform() < 0.5 ? -1.0 : 1.0;
            const double r = sign * std::pow(10.0, -3.0 + 6.0 * uniform());
            bool far = true;
            for (double other : real)
                far = far && std::abs(r - other) > 0.05 * std::max(std::abs(r), std::abs(other));
            if (far)
                real.push_back(r);
        }
        std::vector<std::vector<double>> quadratics;
        for (int i = n_real; i < degree; i += 2) {
            // (x - a)^2 + b^2, roots a +/- i b
            const double u = uniform() - 0.5;
            const double a = u * std::pow(10.0, -2.0 + 4.0 * uniform());
            const double b = std::pow(10.0, -2.0 + 4.0 * uniform());
            quadratics.push_back({a * a + b * b, -2.0 * a, 1.0});
        }
        const std::vector<double> p = poly_from_roots(real, quadratics, 0.5 + uniform());
        CAPTURE(t, degree, n_real);
        const std::vector<double> found = degree == 9 ? scaled_roots<9>(p) : scaled_roots<6>(p);
        check_roots(real, found, 1e-8);
    }
}
