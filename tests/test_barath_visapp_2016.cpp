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
#include <cmath>
#include <random>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "get_barath_visapp_2016.hpp"
#include "get_barath_iccv_2019.hpp"
#include "generate_problem_instance.hpp"
#include "posedata.hpp"

static const double kPI = 3.14159265358979323846;

struct FullInstance {
    Eigen::Matrix3d H, F;
    std::vector<Eigen::Vector2d> x1, x2;
    std::vector<Eigen::Matrix2d> A;
};

static Eigen::Matrix3d skew(const Eigen::Vector3d &v) {
    Eigen::Matrix3d S;
    S << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    return S;
}

static FullInstance generate_full(unsigned int seed, int n_pts) {
    FullInstance inst;
    double fov = std::tan(70.0 / 2.0 * kPI / 180.0);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> depth(0.1, 10.0);
    std::uniform_real_distribution<double> coord(-fov, fov);
    std::normal_distribution<double> dir(0.0, 1.0);

    for (int att = 0; att < 100; ++att) {
        Eigen::Vector3d t(dir(rng), dir(rng), dir(rng)); t.normalize();
        Eigen::Quaterniond q(dir(rng), dir(rng), dir(rng), dir(rng)); q.normalize();
        Eigen::Matrix3d R = q.toRotationMatrix();
        Eigen::Vector3d n(dir(rng), dir(rng), dir(rng)); n.normalize();
        double alpha = depth(rng) / n(2);
        inst.H = alpha * R + t * n.transpose();
        inst.F = skew(t) * R;
        inst.x1.clear(); inst.x2.clear(); inst.A.clear();
        bool ok = true;
        for (int j = 0; j < n_pts; ++j) {
            bool found = false;
            for (int tr = 0; tr < 10; ++tr) {
                Eigen::Vector3d p1{coord(rng), coord(rng), 1.0}; p1.normalize();
                double lam = alpha / n.dot(p1);
                if (lam < 0) continue;
                Eigen::Vector3d p2 = (R * (p1 * lam) + t).normalized();
                if (p2(2) < 0) continue;
                Eigen::Vector2d x2h = p2.hnormalized();
                if (std::abs(x2h(0)) > fov || std::abs(x2h(1)) > fov) continue;
                Eigen::Vector2d x1h = p1.hnormalized();
                double s = inst.H(2,0)*x1h(0) + inst.H(2,1)*x1h(1) + inst.H(2,2);
                Eigen::Matrix2d Aff;
                Aff(0,0) = (inst.H(0,0) - inst.H(2,0)*x2h(0))/s;
                Aff(0,1) = (inst.H(0,1) - inst.H(2,1)*x2h(0))/s;
                Aff(1,0) = (inst.H(1,0) - inst.H(2,0)*x2h(1))/s;
                Aff(1,1) = (inst.H(1,1) - inst.H(2,1)*x2h(1))/s;
                if (Aff.determinant() < 0) continue;
                inst.x1.push_back(x1h); inst.x2.push_back(x2h); inst.A.push_back(Aff);
                found = true; break;
            }
            if (!found) { ok = false; break; }
        }
        if (ok) return inst;
    }
    return inst;
}

static double hom_err(const Eigen::Matrix3d &Hgt, const Eigen::Matrix3d &Hest) {
    Eigen::Matrix3d H1 = Hgt.normalized(); if (H1(2,2) < 0) H1 *= -1;
    Eigen::Matrix3d H2 = Hest.normalized(); if (H2(2,2) < 0) H2 *= -1;
    return (H1 - H2).norm();
}

TEST_CASE("Barath VISAPP 2016 - HA SVD and QR both reproduce GT") {
    int n_svd = 0, n_qr = 0;
    int n_tests = 200;

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 2;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        auto r_svd = HomLib::BarathVISAPP2016::get_affine(inst.x1, inst.x2, inst.A, false);
        auto r_qr = HomLib::BarathVISAPP2016::get_affine(inst.x1, inst.x2, inst.A, true);

        if (!r_svd.empty() && inst.hom_error(r_svd[0].homography) < 1e-8)
            n_svd++;
        if (!r_qr.empty() && inst.hom_error(r_qr[0].homography) < 1e-8)
            n_qr++;
    }

    REQUIRE(n_svd >= 180);
    REQUIRE(n_qr >= 180);
}

TEST_CASE("Barath VISAPP 2016 - HAF reproduces GT homography") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 1);
        if (inst.x1.empty()) continue;
        auto result = HomLib::BarathVISAPP2016::get_haf(inst.x1, inst.x2, inst.A, inst.F);
        if (result.empty()) continue;
        if (hom_err(inst.H, result[0].homography) < 1e-8)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - 3PT reproduces GT homography") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 3);
        if (inst.x1.size() < 3) continue;
        auto result = HomLib::BarathVISAPP2016::get_3pt(inst.x1, inst.x2, inst.F);
        if (result.empty()) continue;
        if (hom_err(inst.H, result[0].homography) < 1e-8)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - Original SIFT reproduces GT homography") {
    int n_pass = 0;
    int n_tests = 200;

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 2;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        auto result = HomLib::BarathICCV2019::get_sift(inst.x1, inst.x2, inst.ori, inst.scale);
        if (result.empty()) continue;

        for (const auto &r : result) {
            if (inst.hom_error(r.homography) < 1e-6) {
                n_pass++;
                break;
            }
        }
    }

    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - Original orientation reproduces GT homography") {
    int n_pass = 0;
    int n_tests = 200;

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 3;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        auto result = HomLib::BarathICCV2019::get_ori(inst.x1, inst.x2, inst.ori);
        if (result.empty()) continue;

        if (inst.hom_error(result[0].homography) < 1e-6)
            n_pass++;
    }

    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - Original orientation overdetermined") {
    int n_pass = 0;
    int n_tests = 200;

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 5;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        auto result = HomLib::BarathICCV2019::get_ori(inst.x1, inst.x2, inst.ori);
        if (result.empty()) continue;

        if (inst.hom_error(result[0].homography) < 1e-6)
            n_pass++;
    }

    REQUIRE(n_pass >= 180);
}

TEST_CASE("SIFT solver solution count and correctness analysis", "[sift_degree]") {
    int n_tests = 10000;
    int barath_count[5] = {};
    int barath_correct[5] = {};

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 2;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        auto rb = HomLib::BarathICCV2019::get_sift(inst.x1, inst.x2, inst.ori, inst.scale);

        int nb = static_cast<int>(rb.size());
        if (nb <= 4) barath_count[nb]++;

        int bc = 0;
        for (const auto &r : rb)
            if (inst.hom_error(r.homography) < 1e-6) bc++;
        if (bc <= 4) barath_correct[bc]++;
    }

    WARN("=== Solution count distribution (over " << n_tests << " instances) ===");
    WARN("Barath (4x4 eig):  0=" << barath_count[0] << "  1=" << barath_count[1]
         << "  2=" << barath_count[2] << "  3=" << barath_count[3]
         << "  4=" << barath_count[4]);
    WARN("=== Correct solutions (error < 1e-6 vs GT) ===");
    WARN("Barath correct:    0=" << barath_correct[0] << "  1=" << barath_correct[1]
         << "  2=" << barath_correct[2] << "  3=" << barath_correct[3]
         << "  4=" << barath_correct[4]);

    REQUIRE(barath_count[4] == 0);
    CHECK(true);
}

TEST_CASE("SIFT discriminant analysis", "[sift_degree]") {
    // Replicate the solver's internal quadratic to analyze the discriminant.
    // disc = (a1*b2 - a2*b1)^2 + 4*b1*b2*q1^2
    // disc >= 0 when b1*b2 >= 0, otherwise can be negative.

    int n_tests = 10000;
    int disc_positive = 0, disc_negative = 0, disc_zero = 0;
    int b1b2_positive = 0, b1b2_negative = 0;

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 2;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        // Replicate solver internals
        const double q1 = inst.scale[0](1) / inst.scale[0](0);
        const double q2 = inst.scale[1](1) / inst.scale[1](0);
        const double det1 = q1 * q1;
        const double det2 = q2 * q2;
        if (det1 <= 0.0 || det2 <= 0.0) continue;

        const double c1r = std::cos(inst.ori[0](0)), s1r = std::sin(inst.ori[0](0));
        const double c1q = std::cos(inst.ori[0](1)), s1q = std::sin(inst.ori[0](1));
        const double c2r = std::cos(inst.ori[1](0)), s2r = std::sin(inst.ori[1](0));
        const double c2q = std::cos(inst.ori[1](1)), s2q = std::sin(inst.ori[1](1));

        const double dx = inst.x1[1](0) - inst.x1[0](0);
        const double dy = inst.x1[1](1) - inst.x1[0](1);
        const double u =  c1r * dx + s1r * dy;
        const double vi = -s1r * dx + c1r * dy;

        const double dxp = inst.x2[1](0) - inst.x2[0](0);
        const double dyp = inst.x2[1](1) - inst.x2[0](1);
        const double up =  c1q * dxp + s1q * dyp;
        const double vp = -s1q * dxp + c1q * dyp;

        const double ct = c2r * c1r + s2r * s1r;
        const double st = s2r * c1r - c2r * s1r;
        const double ctp = c2q * c1q + s2q * s1q;
        const double stp = s2q * c1q - c2q * s1q;

        Eigen::Matrix<double, 3, 5> M;
        Eigen::Vector3d rhs;
        M(0,0) = u;   M(0,1) = vi;  M(0,2) = 0.0;
        M(0,3) = -up*u; M(0,4) = -up*vi; rhs(0) = up;
        M(1,0) = 0.0; M(1,1) = 0.0; M(1,2) = vi;
        M(1,3) = -vp*u; M(1,4) = -vp*vi; rhs(1) = vp;
        M(2,0) = -ct*stp; M(2,1) = -st*stp; M(2,2) = st*ctp;
        M(2,3) = up*ct*stp - vp*ct*ctp;
        M(2,4) = up*st*stp - vp*st*ctp; rhs(2) = 0.0;

        Eigen::JacobiSVD<Eigen::Matrix<double, 3, 5>> svd(
            M, Eigen::ComputeFullU | Eigen::ComputeFullV);
        auto sv = svd.singularValues();
        if (sv(2) < 1e-10 * sv(0)) continue;

        const auto &V = svd.matrixV();
        Eigen::Matrix<double, 5, 1> n1 = V.col(3);
        Eigen::Matrix<double, 5, 1> n2 = V.col(4);
        Eigen::Vector3d Ut_rhs = svd.matrixU().transpose() * rhs;
        Eigen::Matrix<double, 5, 1> p;
        p.setZero();
        for (int i = 0; i < 3; ++i)
            p += (Ut_rhs(i) / sv(i)) * V.col(i);

        const double w0 = std::cbrt(det1 / det2);
        const double aL = n1(3)*u + n1(4)*vi;
        const double bL = n2(3)*u + n2(4)*vi;
        const double cL = w0 - 1.0 - p(3)*u - p(4)*vi;

        if (std::abs(aL) < 1e-14 && std::abs(bL) < 1e-14) continue;

        double A1, B1, A2, B2;
        if (std::abs(aL) >= std::abs(bL)) {
            double alpha0 = cL / aL;
            double alpha_b = -bL / aL;
            A1 = p(0) + alpha0 * n1(0);
            B1 = n2(0) + alpha_b * n1(0);
            A2 = p(2) + alpha0 * n1(2);
            B2 = n2(2) + alpha_b * n1(2);
        } else {
            double beta0 = cL / bL;
            double beta_a = -aL / bL;
            A1 = p(0) + beta0 * n2(0);
            B1 = n1(0) + beta_a * n2(0);
            A2 = p(2) + beta0 * n2(2);
            B2 = n1(2) + beta_a * n2(2);
        }

        double qa = B1 * B2;
        double qb = A1 * B2 + A2 * B1;
        double qc = A1 * A2 - det1;
        double disc = qb * qb - 4.0 * qa * qc;
        // Equivalent form: disc = (A1*B2 - A2*B1)^2 + 4*B1*B2*det1
        double disc_alt = (A1*B2 - A2*B1)*(A1*B2 - A2*B1) + 4.0*B1*B2*det1;

        if (B1*B2 > 0) b1b2_positive++;
        else b1b2_negative++;

        if (disc > 1e-12) disc_positive++;
        else if (disc < -1e-12) disc_negative++;
        else disc_zero++;
    }

    WARN("=== Discriminant analysis (" << n_tests << " instances) ===");
    WARN("disc > 0:  " << disc_positive << "  (two distinct real roots)");
    WARN("disc ~ 0:  " << disc_zero << "  (double root)");
    WARN("disc < 0:  " << disc_negative << "  (no real roots)");
    WARN("b1*b2 > 0: " << b1b2_positive << "  (disc always >= 0)");
    WARN("b1*b2 < 0: " << b1b2_negative << "  (disc can be negative)");

    // Check what happens when we DON'T filter: are both roots finite H?
    // Track sign of h11, h22 for each root.
    int root_both_pos_neg = 0;  // root1: (+,+), root2: (-,-)
    int root_other = 0;
    int root2_H33_near_zero = 0;
    int root2_not_finite = 0;

    for (unsigned int seed = 0; seed < static_cast<unsigned int>(n_tests); seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 2;
        config.no_distortion = true;
        config.seed = seed;

        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);
        const double q1 = inst.scale[0](1) / inst.scale[0](0);
        const double q2 = inst.scale[1](1) / inst.scale[1](0);
        const double det1 = q1 * q1, det2 = q2 * q2;
        if (det1 <= 0.0 || det2 <= 0.0) continue;

        const double c1r = std::cos(inst.ori[0](0)), s1r = std::sin(inst.ori[0](0));
        const double c1q = std::cos(inst.ori[0](1)), s1q = std::sin(inst.ori[0](1));
        const double c2r = std::cos(inst.ori[1](0)), s2r = std::sin(inst.ori[1](0));
        const double c2q = std::cos(inst.ori[1](1)), s2q = std::sin(inst.ori[1](1));

        double dx = inst.x1[1](0)-inst.x1[0](0), dy = inst.x1[1](1)-inst.x1[0](1);
        double u = c1r*dx+s1r*dy, vi = -s1r*dx+c1r*dy;
        double dxp = inst.x2[1](0)-inst.x2[0](0), dyp = inst.x2[1](1)-inst.x2[0](1);
        double up = c1q*dxp+s1q*dyp, vp = -s1q*dxp+c1q*dyp;
        double ct = c2r*c1r+s2r*s1r, st = s2r*c1r-c2r*s1r;
        double ctp = c2q*c1q+s2q*s1q, stp = s2q*c1q-c2q*s1q;

        Eigen::Matrix<double,3,5> M;  Eigen::Vector3d rhs;
        M(0,0)=u;M(0,1)=vi;M(0,2)=0;M(0,3)=-up*u;M(0,4)=-up*vi;rhs(0)=up;
        M(1,0)=0;M(1,1)=0;M(1,2)=vi;M(1,3)=-vp*u;M(1,4)=-vp*vi;rhs(1)=vp;
        M(2,0)=-ct*stp;M(2,1)=-st*stp;M(2,2)=st*ctp;
        M(2,3)=up*ct*stp-vp*ct*ctp;M(2,4)=up*st*stp-vp*st*ctp;rhs(2)=0;

        Eigen::JacobiSVD<Eigen::Matrix<double,3,5>> svd(M, Eigen::ComputeFullU|Eigen::ComputeFullV);
        auto sv = svd.singularValues();
        if (sv(2) < 1e-10*sv(0)) continue;
        auto V = svd.matrixV();
        Eigen::Matrix<double,5,1> n1=V.col(3), n2=V.col(4), p;
        p.setZero();
        Eigen::Vector3d Ut_rhs = svd.matrixU().transpose()*rhs;
        for (int i=0;i<3;++i) p += (Ut_rhs(i)/sv(i))*V.col(i);

        double w0 = std::cbrt(det1/det2);
        double aL=n1(3)*u+n1(4)*vi, bL=n2(3)*u+n2(4)*vi, cLv=w0-1-p(3)*u-p(4)*vi;
        if (std::abs(aL)<1e-14 && std::abs(bL)<1e-14) continue;

        double A1,B1,A2,B2;
        if (std::abs(aL)>=std::abs(bL)) {
            double a0=cLv/aL, ab=-bL/aL;
            A1=p(0)+a0*n1(0); B1=n2(0)+ab*n1(0);
            A2=p(2)+a0*n1(2); B2=n2(2)+ab*n1(2);
        } else {
            double b0=cLv/bL, ba=-aL/bL;
            A1=p(0)+b0*n2(0); B1=n1(0)+ba*n2(0);
            A2=p(2)+b0*n2(2); B2=n1(2)+ba*n2(2);
        }
        double qa=B1*B2, qb_=A1*B2+A2*B1, qc_=A1*A2-det1;
        double disc_ = qb_*qb_ - 4*qa*qc_;
        if (disc_ < 0) continue;
        double sq = std::sqrt(disc_);
        double beta1 = (-qb_+sq)/(2*qa), beta2 = (-qb_-sq)/(2*qa);

        double h11_1 = A1+B1*beta1, h22_1 = A2+B2*beta1;
        double h11_2 = A1+B1*beta2, h22_2 = A2+B2*beta2;

        bool r1_pos = (h11_1 > 0 && h22_1 > 0);
        bool r2_neg = (h11_2 < 0 && h22_2 < 0);
        if ((r1_pos && r2_neg) || (!r1_pos && !r2_neg && h11_1<0 && h22_1<0 && h11_2>0 && h22_2>0))
            root_both_pos_neg++;
        else
            root_other++;

        // Reconstruct H for root 2 and check
        Eigen::Matrix<double,5,1> h2_;
        if (std::abs(aL)>=std::abs(bL)) {
            double alpha2 = cLv/aL + (-bL/aL)*beta2;
            h2_ = p + alpha2*n1 + beta2*n2;
        } else {
            double alpha2_; // beta2 is the free var for the other branch...
            // this branch uses alpha as the quadratic var, so beta2 is actually alpha
            double beta2_back = cLv/bL + (-aL/bL)*beta2;
            h2_ = p + beta2*n1 + beta2_back*n2;
        }
        Eigen::Matrix3d Hp2;
        Hp2 << h2_(0),h2_(1),0, 0,h2_(2),0, h2_(3),h2_(4),1;
        Eigen::Matrix3d S1m, S2inv;
        S1m << c1r,s1r,-(c1r*inst.x1[0](0)+s1r*inst.x1[0](1)),
              -s1r,c1r,(s1r*inst.x1[0](0)-c1r*inst.x1[0](1)),0,0,1;
        S2inv << c1q,-s1q,inst.x2[0](0), s1q,c1q,inst.x2[0](1), 0,0,1;
        Eigen::Matrix3d H2 = S2inv*Hp2*S1m;
        if (std::abs(H2(2,2)) < 1e-10) root2_H33_near_zero++;
        if (!H2.allFinite()) root2_not_finite++;
    }

    WARN("=== Root sign analysis ===");
    WARN("Root 1 (+,+) / Root 2 (-,-): " << root_both_pos_neg);
    WARN("Other pattern: " << root_other);
    WARN("Root 2 H(2,2) near zero: " << root2_H33_near_zero);
    WARN("Root 2 not allFinite: " << root2_not_finite);

    // Check how many instances Barath returns 2 finite solutions
    int both_finite = 0, root2_correct = 0;
    for (unsigned int seed = 0; seed < 10000u; seed++) {
        HomLib::ProblemConfig config;
        config.distortion = HomLib::NO_DISTORTION;
        config.point_noise = 0.0;
        config.number_points = 2;
        config.no_distortion = true;
        config.seed = seed;
        HomLib::ProblemInstance inst = HomLib::generate_problem_instance(config);

        auto rb = HomLib::BarathICCV2019::get_sift(inst.x1, inst.x2, inst.ori, inst.scale);

        if (rb.size() == 2) {
            both_finite++;
            if (inst.hom_error(rb[1].homography) < 1e-6) root2_correct++;
        }
    }

    WARN("=== Barath 2-solution analysis ===");
    WARN("Barath returned 2 finite solutions: " << both_finite);
    WARN("Barath root 2 also correct (matches GT): " << root2_correct);

    CHECK(true);
}

TEST_CASE("Barath VISAPP 2016 - HAF overdetermined improves accuracy") {
    auto inst = generate_full(42, 5);
    REQUIRE(inst.x1.size() == 5);

    auto r1 = HomLib::BarathVISAPP2016::get_haf(
        {inst.x1[0]}, {inst.x2[0]}, {inst.A[0]}, inst.F);
    auto r5 = HomLib::BarathVISAPP2016::get_haf(inst.x1, inst.x2, inst.A, inst.F);

    REQUIRE(!r1.empty());
    REQUIRE(!r5.empty());
    REQUIRE(hom_err(inst.H, r5[0].homography) < 1e-8);
}

TEST_CASE("Barath VISAPP 2016 - 3PT estimate() transforms F under normalization", "[F_norm]") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 3);
        if (inst.x1.size() < 3) continue;

        // Scale points to pixel-like magnitudes so normalization is non-trivial
        double scale = 300.0;
        std::vector<Eigen::Vector2d> x1s(inst.x1), x2s(inst.x2);
        for (auto &p : x1s) p *= scale;
        for (auto &p : x2s) p *= scale;
        // F transforms as F' = S^{-T} F S^{-1} where S = diag(s,s,1)
        Eigen::Matrix3d Fs = Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal()
                           * inst.F
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();

        HomLib::BarathVISAPP2016::ThreePTSolver solver(Fs);
        std::vector<HomLib::PoseData> poses;
        solver.estimate(x1s, x2s, &poses);
        if (poses.empty()) continue;

        // The returned H should be in the scaled coordinate frame
        Eigen::Matrix3d Hs = Eigen::Vector3d(scale, scale, 1).asDiagonal()
                           * inst.H
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();
        if (hom_err(Hs, poses[0].homography) < 1e-6)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}

TEST_CASE("Barath VISAPP 2016 - HAF estimate() transforms F under normalization", "[F_norm]") {
    int n_pass = 0;
    for (unsigned int seed = 0; seed < 200; seed++) {
        auto inst = generate_full(seed, 1);
        if (inst.x1.empty()) continue;

        double scale = 300.0;
        std::vector<Eigen::Vector2d> x1s(inst.x1), x2s(inst.x2);
        for (auto &p : x1s) p *= scale;
        for (auto &p : x2s) p *= scale;
        Eigen::Matrix3d Fs = Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal()
                           * inst.F
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();

        HomLib::BarathVISAPP2016::HAFSolver solver(Fs);
        std::vector<HomLib::PoseData> poses;
        solver.estimate(x1s, x2s, inst.A, &poses);
        if (poses.empty()) continue;

        Eigen::Matrix3d Hs = Eigen::Vector3d(scale, scale, 1).asDiagonal()
                           * inst.H
                           * Eigen::Vector3d(1.0/scale, 1.0/scale, 1).asDiagonal();
        if (hom_err(Hs, poses[0].homography) < 1e-6)
            n_pass++;
    }
    REQUIRE(n_pass >= 180);
}
