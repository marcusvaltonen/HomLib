#pragma once
// The gamma tables of MAGSAC++ (danini/magsac's magsac.h) for 2 degrees of freedom, in place of
// GC-RANSAC's gamma_values.cpp, whose tables are those of 4. HomLib's residual is the 2-D transfer
// error in image 2 (HomLibMAGSACEstimator::getDegreesOfFreedom()). CMake's fixed copy of magsac.h
// includes this file and reads these tables instead of gamma_values.cpp's, which stays as it is
// (it belongs to the GC-RANSAC submodule): the two that magsac.h reads once its IRLS weight is
// fixed (danini/magsac#49), with the same size and step:
// - complete[i]: the upper incomplete gamma Gamma((n-1)/2, i / precision);
// - lower[i]: the lower incomplete gamma gamma((n+1)/2, i / precision).
// With n = 2 both have closed forms: Gamma(1/2, t) = sqrt(pi) erfc(sqrt(t)) and
// gamma(3/2, t) = sqrt(pi)/2 erf(sqrt(t)) - sqrt(t) exp(-t).

#include <cmath>
#include <cstddef>
#include <vector>

namespace HomLib::magsac_gamma {

constexpr double sqrt_pi = 1.7724538509055160273;

inline double upper_half(double t) { return sqrt_pi * std::erfc(std::sqrt(t)); }
inline double lower_three_halves(double t) {
    return sqrt_pi / 2 * std::erf(std::sqrt(t)) - std::sqrt(t) * std::exp(-t);
}

inline std::vector<double> table(double (*f)(double), std::size_t last, double step) {
    std::vector<double> values(last + 1);
    for (std::size_t i = 0; i <= last; ++i)
        values[i] = f(i * step);
    return values;
}

constexpr unsigned int number = 100000;   // the last index
constexpr double precision = 1e4;

inline const std::vector<double> complete = table(upper_half, number, 1 / precision);
inline const std::vector<double> lower = table(lower_three_halves, number, 1 / precision);

} // namespace HomLib::magsac_gamma
