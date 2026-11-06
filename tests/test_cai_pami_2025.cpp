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
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_cai_pami_2025.hpp"
#include "get_dlt.hpp"
#include "get_guo_prl_2019.hpp"
#include "posedata.hpp"

TEST_CASE("Cai PAMI 2025 - basic") {
    // Ground-truth homography
    Eigen::Matrix3d H_gt;
    H_gt << 1.2, 0.3, 0.1,
           -0.1, 0.9, 0.2,
            0.05, -0.02, 1.0;

    // Source points (4 points in general position)
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.5, -0.3),
        Eigen::Vector2d( 0.4,  0.2),
        Eigen::Vector2d(-0.2,  0.5),
        Eigen::Vector2d( 0.3, -0.4)
    };

    // Project through H_gt to get target points
    std::vector<Eigen::Vector2d> p2;
    for (int i = 0; i < 4; i++) {
        Eigen::Vector3d tmp = H_gt * p1[i].homogeneous();
        p2.push_back(tmp.hnormalized());
    }

    std::vector<HomLib::PoseData> posedata = HomLib::CaiPAMI2025::get(p1, p2);

    REQUIRE(posedata.size() == 1);

    // No distortion
    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.0).margin(1e-12));

    // The recovered H should match H_gt up to scale
    Eigen::Matrix3d H_est = posedata[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_ref = H_gt / H_gt(2, 2);
    REQUIRE(H_est.isApprox(H_ref, 1e-10));
}

TEST_CASE("Cai PAMI 2025 - identity homography") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-1.0,  0.0),
        Eigen::Vector2d( 1.0,  0.0),
        Eigen::Vector2d( 0.0,  1.0),
        Eigen::Vector2d( 0.0, -1.0)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::CaiPAMI2025::get(p1, p1);

    REQUIRE(posedata.size() == 1);

    Eigen::Matrix3d H_est = posedata[0].homography;
    H_est /= H_est(2, 2);
    REQUIRE(H_est.isApprox(Eigen::Matrix3d::Identity(), 1e-10));
}

TEST_CASE("Cai PAMI 2025 - degenerate (collinear points)") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(0.0, 0.0),
        Eigen::Vector2d(1.0, 0.0),
        Eigen::Vector2d(2.0, 0.0),
        Eigen::Vector2d(3.0, 0.0)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::CaiPAMI2025::get(p1, p1);

    REQUIRE(posedata.size() == 0);
}

TEST_CASE("Cai PAMI 2025 - perspective homography") {
    Eigen::Matrix3d H_gt;
    H_gt << 0.8, -0.2, 0.5,
            0.1,  1.1, -0.3,
           -0.1,  0.05, 1.0;

    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.4, -0.6),
        Eigen::Vector2d( 0.5,  0.3),
        Eigen::Vector2d(-0.3,  0.4),
        Eigen::Vector2d( 0.6, -0.2)
    };

    std::vector<Eigen::Vector2d> p2;
    for (int i = 0; i < 4; i++) {
        Eigen::Vector3d tmp = H_gt * p1[i].homogeneous();
        p2.push_back(tmp.hnormalized());
    }

    std::vector<HomLib::PoseData> posedata = HomLib::CaiPAMI2025::get(p1, p2);

    REQUIRE(posedata.size() == 1);

    Eigen::Matrix3d H_est = posedata[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_ref = H_gt / H_gt(2, 2);
    REQUIRE(H_est.isApprox(H_ref, 1e-10));
}

// ---- DLT tests ----

TEST_CASE("DLT - basic") {
    Eigen::Matrix3d H_gt;
    H_gt << 1.2, 0.3, 0.1,
           -0.1, 0.9, 0.2,
            0.05, -0.02, 1.0;

    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.5, -0.3),
        Eigen::Vector2d( 0.4,  0.2),
        Eigen::Vector2d(-0.2,  0.5),
        Eigen::Vector2d( 0.3, -0.4)
    };

    std::vector<Eigen::Vector2d> p2;
    for (int i = 0; i < 4; i++) {
        Eigen::Vector3d tmp = H_gt * p1[i].homogeneous();
        p2.push_back(tmp.hnormalized());
    }

    std::vector<HomLib::PoseData> posedata = HomLib::DLT::get(p1, p2);

    REQUIRE(posedata.size() == 1);

    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.0).margin(1e-12));

    Eigen::Matrix3d H_est = posedata[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_ref = H_gt / H_gt(2, 2);
    REQUIRE(H_est.isApprox(H_ref, 1e-10));
}

TEST_CASE("DLT - degenerate (collinear points)") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(0.0, 0.0),
        Eigen::Vector2d(1.0, 0.0),
        Eigen::Vector2d(2.0, 0.0),
        Eigen::Vector2d(3.0, 0.0)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::DLT::get(p1, p1);

    REQUIRE(posedata.size() == 0);
}

TEST_CASE("DLT QR - basic") {
    Eigen::Matrix3d H_gt;
    H_gt << 1.2, 0.3, 0.1,
           -0.1, 0.9, 0.2,
            0.05, -0.02, 1.0;

    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.5, -0.3),
        Eigen::Vector2d( 0.4,  0.2),
        Eigen::Vector2d(-0.2,  0.5),
        Eigen::Vector2d( 0.3, -0.4)
    };

    std::vector<Eigen::Vector2d> p2;
    for (int i = 0; i < 4; i++) {
        Eigen::Vector3d tmp = H_gt * p1[i].homogeneous();
        p2.push_back(tmp.hnormalized());
    }

    std::vector<HomLib::PoseData> posedata = HomLib::DLT::get(p1, p2, true);

    REQUIRE(posedata.size() == 1);

    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.0).margin(1e-12));

    Eigen::Matrix3d H_est = posedata[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_ref = H_gt / H_gt(2, 2);
    REQUIRE(H_est.isApprox(H_ref, 1e-10));
}

// ---- Guo PRL 2019 tests ----

TEST_CASE("Guo PRL 2019 - basic") {
    Eigen::Matrix3d H_gt;
    H_gt << 1.2, 0.3, 0.1,
           -0.1, 0.9, 0.2,
            0.05, -0.02, 1.0;

    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-0.5, -0.3),
        Eigen::Vector2d( 0.4,  0.2),
        Eigen::Vector2d(-0.2,  0.5),
        Eigen::Vector2d( 0.3, -0.4)
    };

    std::vector<Eigen::Vector2d> p2;
    for (int i = 0; i < 4; i++) {
        Eigen::Vector3d tmp = H_gt * p1[i].homogeneous();
        p2.push_back(tmp.hnormalized());
    }

    std::vector<HomLib::PoseData> posedata = HomLib::GuoPRL2019::get(p1, p2);

    REQUIRE(posedata.size() == 1);

    REQUIRE(posedata[0].distortion_parameter == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(posedata[0].distortion_parameter2 == Catch::Approx(0.0).margin(1e-12));

    Eigen::Matrix3d H_est = posedata[0].homography;
    H_est /= H_est(2, 2);
    Eigen::Matrix3d H_ref = H_gt / H_gt(2, 2);
    REQUIRE(H_est.isApprox(H_ref, 1e-10));
}

TEST_CASE("Guo PRL 2019 - identity homography") {
    std::vector<Eigen::Vector2d> p1 = {
        Eigen::Vector2d(-1.0,  0.0),
        Eigen::Vector2d( 1.0,  0.0),
        Eigen::Vector2d( 0.0,  1.0),
        Eigen::Vector2d( 0.0, -1.0)
    };

    std::vector<HomLib::PoseData> posedata = HomLib::GuoPRL2019::get(p1, p1);

    REQUIRE(posedata.size() == 1);

    Eigen::Matrix3d H_est = posedata[0].homography;
    H_est /= H_est(2, 2);
    REQUIRE(H_est.isApprox(Eigen::Matrix3d::Identity(), 1e-10));
}
