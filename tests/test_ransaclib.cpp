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

// Tests of RansacLib's NumRequiredIterations, through the fixed copy of RansacLib/utils.h that
// CMakeLists.txt writes (tsattler/RansacLib#24).

#include <RansacLib/utils.h>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <vector>

using ransac_lib::utils::NumRequiredIterations;

namespace {
// RansacLib's default success probability is 0.9999 (ransac.h), so the probability of missing
// the best model is 1e-4. Sample size 4: a homography from four points.
constexpr double kProbMissing = 1e-4;
constexpr int kSampleSize = 4;
constexpr uint32_t kMinIterations = 100;
constexpr uint32_t kMaxIterations = 10000;
}  // namespace

TEST_CASE("RansacLib - an iteration count above 2^32 gives the maximum") {
    // About 2^32 + 381 iterations are needed: converting the count to uint32_t before clamping
    // it gives 381 on x86-64. (volatile: computed at run time, as in RANSAC, not folded by the
    // compiler, which may give anything for an undefined conversion.)
    const volatile double inlier_ratio = 0.006805012733;
    REQUIRE(NumRequiredIterations(inlier_ratio, kProbMissing, kSampleSize, kMinIterations,
                                  kMaxIterations) == kMaxIterations);
    REQUIRE(NumRequiredIterations(std::vector<double>{inlier_ratio}, kProbMissing,
                                  std::vector<int>{kSampleSize}, kMinIterations,
                                  kMaxIterations) == kMaxIterations);
}

TEST_CASE("RansacLib - a NaN iteration count gives the maximum") {
    // Fails with RansacLib's own utils.h on any platform: converting NaN to uint32_t gives 0 on
    // x86-64 and ARM64, which becomes min_iterations. So this also checks that the fixed copy is
    // the one included.
    const volatile double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(NumRequiredIterations(nan, kProbMissing, kSampleSize, kMinIterations, kMaxIterations) ==
            kMaxIterations);
    REQUIRE(NumRequiredIterations(std::vector<double>{nan}, kProbMissing, std::vector<int>{kSampleSize},
                                  kMinIterations, kMaxIterations) == kMaxIterations);
}

TEST_CASE("RansacLib - small iteration counts are unchanged") {
    // ceil(log(1e-4) / log(1 - 0.2^4) + 0.5) = 5753, between the minimum and the maximum
    REQUIRE(NumRequiredIterations(0.2, kProbMissing, kSampleSize, kMinIterations, kMaxIterations) == 5753);
    // 0.9^4 = 0.6561: 10 iterations, raised to the minimum
    REQUIRE(NumRequiredIterations(0.9, kProbMissing, kSampleSize, kMinIterations, kMaxIterations) ==
            kMinIterations);
}
