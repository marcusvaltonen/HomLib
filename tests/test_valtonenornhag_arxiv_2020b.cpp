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
#include <catch2/catch.hpp>
#include "get_valtonenornhag_arxiv_2020b.hpp"
#include "posedata.hpp"

TEST_CASE("Valtonen Ornhag Arxiv 2020 B - fHf") {
    Eigen::MatrixXd p1(2, 2);
    Eigen::MatrixXd p2(2, 2);
    Eigen::Matrix3d R1, R2;
    p1 << 2.003107199098924, -15.634084933471335,
         -0.017087350257598,  -7.041596829586987,
    p2 << 0.395688457559412,  -0.012777594199286,
          2.097270018093999,   0.988175585551782;
    R1 << 0.854451801803156, 0.080889542675225,  0.513194895026376,
          0.251645807638113, 0.799754574299643, -0.545038538440134,
         -0.454517882919365, 0.594852505057001,  0.662996222714662;

    R2 << 0.243935353955667, -0.895887070857591, -0.371324520279403,
          0.945623784801441,  0.134783533871648,  0.296022054271079,
         -0.215153900053711, -0.423343542843488,  0.880050591741408;

    std::vector<HomLib::PoseData> posedata = HomLib::ValtonenOrnhagArxiv2020B::get_fHf(p1, p2, R1, R2);

    double tol = 1e-12;

    // Test size
    REQUIRE(posedata.size() == 4);

    // Test distortion parameters
    REQUIRE(posedata[0].focal_length == Approx(-0.169636953030949).margin(tol));
    REQUIRE(posedata[1].focal_length == Approx(0.341654151554244).margin(tol));
    REQUIRE(posedata[2].focal_length == Approx(0.653592155944414).margin(tol));
    REQUIRE(posedata[3].focal_length == Approx(3.395609550168773).margin(tol));

    // Test homographies
    tol = 1e-9;
    Eigen::Matrix3d expected;

    // Test first putative homography
    expected << 0.000942286023272176, -0.00406723298360957,  -0.0139080351631117,
                  -0.119952192496138,     0.26057885034196,  -0.0404515590157271,
                   0.236884094847984,   -0.539799164416777,  -0.0975820265293806;
    REQUIRE(posedata[0].homography.isApprox(expected, tol));

    // Test second putative homography
    expected << 0.0119880554013824,   0.145480867162921, -0.0215274775447118,
                 0.165941018970929,   -1.27273026931389,   -0.35414523762168,
                0.0803816918666426,   -2.08713777622907,  -0.196676799872108;
    REQUIRE(posedata[1].homography.isApprox(expected, tol));

    // Test third putative homography
    expected << -0.00407848497852057,   0.0112415881437603,   -0.208614339370427,
                    0.23741442102905,     -3.2260628888446,    -1.68073103477477,
                  0.0248475548331138,     -2.6362035596848,   -0.643168820965163;
    REQUIRE(posedata[2].homography.isApprox(expected, tol));

    // Test fourth putative homography
    expected << -0.104788674142645, -0.768146499303879,  -7.71830881514091,
                  1.25789148717968,  -16.5186232373972,  -44.7543319830629,
                -0.600283285536537,  -8.81686959619256,  -18.9515536030011;
    REQUIRE(posedata[3].homography.isApprox(expected, tol));
}

TEST_CASE("Valtonen Ornhag Arxiv 2020 B - frHfr") {
    Eigen::MatrixXd p1(2, 3);
    Eigen::MatrixXd p2(2, 3);
    Eigen::Matrix3d R1, R2;

    p1 << -1.146281331283839,   1.050109134098990,   1.065996624259908,
          -0.879951300627967,   0.620743795713172,   0.541580087112080;
    p2 <<  0.663628650450811,   1.333268512835822,   1.318951998842419,
           1.241359691717976,   0.068745345721370,   0.016786262835316;
    R1 << -0.320761154096478,   0.935110133718446,   0.150603186685294,
          -0.808554515336552,  -0.353152142517100,   0.470662469254195,
           0.493307082608358,   0.029199350209406,   0.869365009760445;
    R2 <<  0.420952545706761,  -0.744893719609945,  -0.517621773835547,
           0.652124812245433,  -0.148125936460580,   0.743499788972085,
          -0.630501533318403,  -0.650532130976893,   0.423409687005158;


    HomLib::PoseData posedata = HomLib::ValtonenOrnhagArxiv2020B::get_frHfr(p1, p2, R1, R2);

    double tol = 1e-12;

    // Test distortion parameter and focal length
    REQUIRE(posedata.distortion_parameter == Approx(-0.4505415985025416).margin(tol));
    REQUIRE(posedata.focal_length == Approx(5.756798219567954).margin(tol));

    // Test homography
    tol = 1e-7;
    Eigen::Matrix3d expected;
    expected <<  -5.13410573,  -2.26448940,  -9.50684998,
                -12.94693678,   3.53534100,  33.04037699,
                 -2.42428225,   1.80519293,  -0.12730716;
    REQUIRE(posedata.homography.isApprox(expected, tol));
}
