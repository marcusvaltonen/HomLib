#ifndef TEST_GENERATE_HELPERS_HPP_
#define TEST_GENERATE_HELPERS_HPP_

#include <Eigen/Dense>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include "radial.hpp"

// Generate noiseless one-sided left correspondences: y is distorted, x is not.
inline void generate_one_sided_left(
    const Eigen::Matrix3d& H, double k2, int n,
    std::vector<Eigen::Vector2d>* x, std::vector<Eigen::Vector2d>* y,
    unsigned seed = 42)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> coord(-0.8, 0.8);
    x->resize(n);
    y->resize(n);
    for (int i = 0; i < n; i++) {
        (*x)[i] = Eigen::Vector2d(coord(rng), coord(rng));
        Eigen::Vector3d tmp = H * (*x)[i].homogeneous();
        Eigen::Vector2d yu = tmp.hnormalized();
        Eigen::Vector2d yd;
        REQUIRE(HomLib::try_radialdistort(yu, k2, yd));
        (*y)[i] = yd;
    }
}

// Generate noiseless two-sided equal correspondences: both x and y are distorted with the same k.
inline void generate_two_sided_equal(
    const Eigen::Matrix3d& H, double k, int n,
    std::vector<Eigen::Vector2d>* x, std::vector<Eigen::Vector2d>* y,
    unsigned seed = 42)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> coord(-0.8, 0.8);
    x->resize(n);
    y->resize(n);
    for (int i = 0; i < n; i++) {
        Eigen::Vector2d xu(coord(rng), coord(rng));
        Eigen::Vector2d xd, yd;
        REQUIRE(HomLib::try_radialdistort(xu, k, xd));
        (*x)[i] = xd;
        Eigen::Vector3d tmp = H * xu.homogeneous();
        Eigen::Vector2d yu = tmp.hnormalized();
        REQUIRE(HomLib::try_radialdistort(yu, k, yd));
        (*y)[i] = yd;
    }
}

#endif  // TEST_GENERATE_HELPERS_HPP_
