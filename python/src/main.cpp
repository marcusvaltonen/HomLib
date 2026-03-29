#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <Eigen/Dense>

#include <RansacLib/ransac.h>

#include <posedata.hpp>
#include <get_fitzgibbon_cvpr_2001.hpp>
#include <get_kukelova_cvpr_2015.hpp>
#include <get_nakano_icpr_2025.hpp>
#include <get_wadenback_3dv_2026.hpp>
#include <ransac_estimator.h>


namespace py = pybind11;
using namespace pybind11::literals;

void preprocess(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    std::vector<Eigen::Vector2d> *x,
    std::vector<Eigen::Vector2d> *y
) {
    if (x_.cols() != y_.cols()) {
        throw std::invalid_argument("x and y should be of equal size.");
    }
    for (Eigen::Index i=0; i < x_.cols(); i++) {
        x->push_back(x_.col(i));
        y->push_back(y_.col(i));
    }
}


std::vector<HomLib::PoseData> estimate_fitzgibbon_cvpr_2001_one_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::FitzgibbonCVPR2001::get_single(x, y);
    return poses;
} 

std::vector<HomLib::PoseData> estimate_fitzgibbon_cvpr_2001_two_sided_equal_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::FitzgibbonCVPR2001::get(x, y);
    return poses;
}

std::vector<HomLib::PoseData> estimate_kukelova_cvpr_2015_two_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    bool dist_equal
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::KukelovaCVPR2015::get(x, y, dist_equal);
    return poses;
}

std::vector<HomLib::PoseData> estimate_kukelova_cvpr_2015_two_sided_6pt_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    bool dist_equal
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::KukelovaCVPR2015::get_6pt(x, y, dist_equal);
    return poses;
}

std::vector<HomLib::PoseData> estimate_nakano_icpr_2025_one_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    bool extra_check
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::NakanoICPR2025::get(x, y, extra_check);
    return poses;
}

std::vector<HomLib::PoseData> estimate_wadenback_3dv_2026_one_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    bool extra_check
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::Wadenback3DV2026::get_one_sided(x, y, extra_check);
    return poses;
}

std::vector<HomLib::PoseData> estimate_wadenback_3dv_2026_two_sided_equal_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    bool extra_check
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::Wadenback3DV2026::get_double_sided_equal(x, y, extra_check);
    return poses;
}

std::vector<HomLib::PoseData> estimate_wadenback_3dv_2026_two_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    bool extra_check
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    std::vector<HomLib::PoseData> poses = HomLib::Wadenback3DV2026::get_double_sided(x, y, extra_check);
    return poses;
}

template <typename Estimator> std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_wrapper(
    Estimator* estimator,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);
    
    ransac_lib::RansacStatistics ransac_stats;
    int inliers = 0;
    HomLib::PoseData best_model;
    
    HomLib::RansacEstimator<Estimator> solver(x, y, *estimator);
    ransac_lib::LocallyOptimizedMSAC<
        HomLib::PoseData,
        std::vector<HomLib::PoseData>,
        HomLib::RansacEstimator<Estimator>> lomsac;
    inliers = lomsac.EstimateModel(options, solver, &best_model, &ransac_stats);
    (void)inliers;  // Suppress warnings

    return std::make_tuple(std::move(best_model), std::move(ransac_stats));
}


std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_fitzgibbon_cvpr_2001_one_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::FitzgibbonCVPR2001::SolverSingleSided estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_fitzgibbon_cvpr_2001_two_sided_equal_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::FitzgibbonCVPR2001::SolverTwoSidedEqual estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_kukelova_cvpr_2015_two_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::KukelovaCVPR2015::SolverTwoSided estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_kukelova_cvpr_2015_two_sided_equal_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::KukelovaCVPR2015::SolverTwoSidedEqual estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_kukelova_cvpr_2015_two_sided_equal_6pt_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::KukelovaCVPR2015::SolverTwoSidedEqual6Pt estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_kukelova_cvpr_2015_two_sided_6pt_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::KukelovaCVPR2015::SolverTwoSided6Pt estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_nakano_icpr_2025_one_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::NakanoICPR2025::SolverSingleSided estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_wadenback_3dv_2026_one_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::Wadenback3DV2026::SolverSingleSided estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_wadenback_3dv_2026_two_sided_equal_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::Wadenback3DV2026::SolverTwoSidedEqual estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}

std::tuple<HomLib::PoseData, ransac_lib::RansacStatistics> lomsac_wadenback_3dv_2026_two_sided_wrapper(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const ransac_lib::LORansacOptions &options
) {
    HomLib::Wadenback3DV2026::SolverTwoSided estimator;
    auto output = lomsac_wrapper(&estimator, x_, y_, options);
    return output;
}


PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Wrappers for C++ core algorithms of HomLib.
        -----------------------

        .. currentmodule:: _core

        .. autosummary::
           :toctree: _generate
           
            PoseData

            estimate_fitzgibbon_cvpr_2001_one_sided
            estimate_fitzgibbon_cvpr_2001_two_sided_equal
            estimate_kukelova_cvpr_2015_two_sided
            estimate_kukelova_cvpr_2015_two_sided_6pt
            estimate_nakano_icpr_2025_one_sided
            estimate_wadenback_3dv_2026_one_sided
            estimate_wadenback_3dv_2026_two_sided_equal
            estimate_wadenback_3dv_2026_two_sided

            lomsac_fitzgibbon_cvpr_2001_one_sided
            lomsac_fitzgibbon_cvpr_2001_two_sided_equal
            lomsac_kukelova_cvpr_2015_two_sided
            lomsac_kukelova_cvpr_2015_two_sided_equal
            lomsac_kukelova_cvpr_2015_two_sided_equal_6pt
            lomsac_kukelova_cvpr_2015_two_sided_6pt
            lomsac_nakano_icpr_2025_one_sided
            lomsac_wadenback_3dv_2026_one_sided
            lomsac_wadenback_3dv_2026_two_sided_equal
            lomsac_wadenback_3dv_2026_two_sided
            
            LORansacOptions
            RansacStatistics
           
    )pbdoc";
    
    py::class_<HomLib::PoseData>(m, "PoseData")
        .def(
            py::init<Eigen::Matrix3d, double, double, double>(),
            "Constructor for PoseData.",
            "homography"_a,
            "focal_length"_a,
            "distortion_parameter"_a,
            "distortion_parameter2"_a
        )
        .def_readwrite("homography", &HomLib::PoseData::homography)
        .def_readwrite("focal_length", &HomLib::PoseData::focal_length)
        .def_readwrite("distortion_parameter", &HomLib::PoseData::distortion_parameter)
        .def_readwrite("distortion_parameter2", &HomLib::PoseData::distortion_parameter2)
        .def("__repr__",
            [](const HomLib::PoseData &p) {
                return "PoseData("
                    "H=[3x3 np.array], "
                    "focal_length=" + std::to_string(p.focal_length) + ", "
                    "distortion_parameter=" + std::to_string(p.distortion_parameter) + ", "
                    "distortion_parameter2=" + std::to_string(p.distortion_parameter2) +
                    ")";
            }
        );
    
    py::class_<ransac_lib::LORansacOptions>(m, "LORansacOptions")
        .def(
            py::init<>()
        )
        .def_readwrite("min_num_iterations", &ransac_lib::LORansacOptions::min_num_iterations_)
        .def_readwrite("max_num_iterations", &ransac_lib::LORansacOptions::max_num_iterations_)
        .def_readwrite("success_probability", &ransac_lib::LORansacOptions::success_probability_)
        .def_readwrite("squared_inlier_threshold", &ransac_lib::LORansacOptions::squared_inlier_threshold_)
        .def_readwrite("random_seed", &ransac_lib::LORansacOptions::random_seed_)
        .def_readwrite("num_lo_steps", &ransac_lib::LORansacOptions::num_lo_steps_)
        .def_readwrite("threshold_multiplier", &ransac_lib::LORansacOptions::threshold_multiplier_)
        .def_readwrite("num_lsq_iterations", &ransac_lib::LORansacOptions::num_lsq_iterations_)
        .def_readwrite("min_sample_multiplicator", &ransac_lib::LORansacOptions::min_sample_multiplicator_)
        .def_readwrite("non_min_sample_multiplier", &ransac_lib::LORansacOptions::non_min_sample_multiplier_)
        .def_readwrite("lo_starting_iterations", &ransac_lib::LORansacOptions::lo_starting_iterations_)
        .def_readwrite("final_least_squares", &ransac_lib::LORansacOptions::final_least_squares_)
        .def("__repr__",
            [](const ransac_lib::LORansacOptions &opt) {
                return "LORansacOptions("
                    "min_num_iterations=" + std::to_string(opt.min_num_iterations_) + ", "
                    "max_num_iterations=" + std::to_string(opt.max_num_iterations_) + ", "
                    "success_probability=" + std::to_string(opt.success_probability_) + ", "
                    "squared_inlier_threshold=" + std::to_string(opt.squared_inlier_threshold_) + ", "
                    "random_seed=" + std::to_string(opt.random_seed_) + ", "
                    "num_lo_steps=" + std::to_string(opt.num_lo_steps_) + ", "
                    "threshold_multiplier=" + std::to_string(opt.threshold_multiplier_) + ", "
                    "num_lsq_iterations=" + std::to_string(opt.num_lsq_iterations_) + ", "
                    "min_sample_multiplicator=" + std::to_string(opt.min_sample_multiplicator_) + ", "
                    "non_min_sample_multiplier=" + std::to_string(opt.non_min_sample_multiplier_) + ", "
                    "lo_starting_iterations=" + std::to_string(opt.lo_starting_iterations_) + ", "
                    "final_least_squares=" + std::to_string(opt.final_least_squares_) +
                    ")";
            }
        );
    
    py::class_<ransac_lib::RansacStatistics>(m, "RansacStatistics")
        .def(
            py::init<>()
        )
        .def_readwrite("num_iterations", &ransac_lib::RansacStatistics::num_iterations)
        .def_readwrite("best_num_inliers", &ransac_lib::RansacStatistics::best_num_inliers)
        .def_readwrite("best_model_score", &ransac_lib::RansacStatistics::best_model_score)
        .def_readwrite("inlier_ratio", &ransac_lib::RansacStatistics::inlier_ratio)
        .def_readwrite("inlier_indices", &ransac_lib::RansacStatistics::inlier_indices)
        .def_readwrite("number_lo_iterations", &ransac_lib::RansacStatistics::number_lo_iterations)
        .def("__repr__",
            [](const ransac_lib::RansacStatistics &s) {
                return "RansacStatistics("
                    "num_iterations=" + std::to_string(s.num_iterations) + ", "
                    "best_num_inliers=" + std::to_string(s.best_num_inliers) + ", "
                    "best_model_score=" + std::to_string(s.best_model_score) + ", "
                    "inlier_ratio=" + std::to_string(s.inlier_ratio) + ", "
                    "inlier_indices=[int list of length " + std::to_string(s.inlier_indices.size()) + "], "
                    "number_lo_iterations=" + std::to_string(s.number_lo_iterations) +
                    ")";
            }
        );
    
    /*
     *  Wrappers for minimal solvers
     */
    m.def(
        "estimate_fitzgibbon_cvpr_2001_one_sided",
        &estimate_fitzgibbon_cvpr_2001_one_sided_wrapper,
        R"pbdoc(
            Solver from [1] as modified by [2].
            
            [1] A. W. Fitzgibbon, "Simultaneous linear estimation of multiple view geometry and lens distortion," In
                Computer Vision and Pattern Recognition (CVPR), 2001
            [2]
        )pbdoc",
        "x"_a,
        "y"_a
    );
    m.def(
        "estimate_fitzgibbon_cvpr_2001_two_sided_equal",
        &estimate_fitzgibbon_cvpr_2001_two_sided_equal_wrapper,
        R"pbdoc(
            Solver from [1].
            
            [1] A. W. Fitzgibbon, "Simultaneous linear estimation of multiple view geometry and lens distortion," In
                Computer Vision and Pattern Recognition (CVPR), 2001
        )pbdoc",
        "x"_a,
        "y"_a
    );
    m.def(
        "estimate_kukelova_cvpr_2015_two_sided",
        &estimate_kukelova_cvpr_2015_two_sided_wrapper,
        R"pbdoc(
            Solver from [1].
            
            [1] Z. Kukelova, J. Heller, M. Bujnak and T. Pajdla, "Radial distortion homography," In Computer Vision
                and Pattern Recognition (CVPR), 2015
        )pbdoc",
        "x"_a,
        "y"_a,
        "dist_equal"_a
    );
    m.def(
        "estimate_kukelova_cvpr_2015_two_sided_6pt",
        &estimate_kukelova_cvpr_2015_two_sided_6pt_wrapper,
        R"pbdoc(
            Solver from [1].
            
            [1] Z. Kukelova, J. Heller, M. Bujnak and T. Pajdla, "Radial distortion homography," In Computer Vision
                and Pattern Recognition (CVPR), 2015
        )pbdoc",
        "x"_a,
        "y"_a,
        "dist_equal"_a
    );
    m.def(
        "estimate_nakano_icpr_2025_one_sided",
        &estimate_nakano_icpr_2025_one_sided_wrapper,
        R"pbdoc(
            Solver from [1]
            
            [1] Gaku Nakano. "Inverse DLT Method for One-Sided Radial Distortion Homography", In
                International Conference on Pattern Recognition (ICPR), 2024.
        )pbdoc",
        "x"_a,
        "y"_a,
        "extra_check"_a
    );
    m.def(
        "estimate_wadenback_3dv_2026_one_sided",
        &estimate_wadenback_3dv_2026_one_sided_wrapper,
        R"pbdoc(
            Solver from [1]
            
            [1] TODO
        )pbdoc",
        "x"_a,
        "y"_a,
        "extra_check"_a
    );
    m.def(
        "estimate_wadenback_3dv_2026_two_sided_equal",
        &estimate_wadenback_3dv_2026_two_sided_equal_wrapper,
        R"pbdoc(
            Solver from [1]
            
            [1] TODO
        )pbdoc",
        "x"_a,
        "y"_a,
        "extra_check"_a
    );
    m.def(
        "estimate_wadenback_3dv_2026_two_sided",
        &estimate_wadenback_3dv_2026_two_sided_wrapper,
        R"pbdoc(
            Solver from [1]
            
            [1] TODO
        )pbdoc",
        "x"_a,
        "y"_a,
        "extra_check"_a
    );
    
    /*
     *  LOMSAC wrappers
     */
    
    m.def(
        "lomsac_fitzgibbon_cvpr_2001_one_sided",
        &lomsac_fitzgibbon_cvpr_2001_one_sided_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] A. W. Fitzgibbon, "Simultaneous linear estimation of multiple view geometry and lens distortion," In
                Conference on Computer Vision and Pattern Recognition (CVPR), 2001
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_fitzgibbon_cvpr_2001_two_sided_equal",
        &lomsac_fitzgibbon_cvpr_2001_two_sided_equal_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] A. W. Fitzgibbon, "Simultaneous linear estimation of multiple view geometry and lens distortion," In
                Conference on Computer Vision and Pattern Recognition (CVPR), 2001
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_kukelova_cvpr_2015_two_sided",
        &lomsac_kukelova_cvpr_2015_two_sided_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] Z. Kukelova, J. Heller, M. Bujnak and T. Pajdla, "Radial distortion homography," In Computer Vision
                and Pattern Recognition (CVPR), 2015
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_kukelova_cvpr_2015_two_sided_equal",
        &lomsac_kukelova_cvpr_2015_two_sided_equal_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] Z. Kukelova, J. Heller, M. Bujnak and T. Pajdla, "Radial distortion homography," In Computer Vision
                and Pattern Recognition (CVPR), 2015
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_kukelova_cvpr_2015_two_sided_equal_6pt",
        &lomsac_kukelova_cvpr_2015_two_sided_equal_6pt_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] Z. Kukelova, J. Heller, M. Bujnak and T. Pajdla, "Radial distortion homography," In Computer Vision
                and Pattern Recognition (CVPR), 2015
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_kukelova_cvpr_2015_two_sided_6pt",
        &lomsac_kukelova_cvpr_2015_two_sided_6pt_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] Z. Kukelova, J. Heller, M. Bujnak and T. Pajdla, "Radial distortion homography," In Computer Vision
                and Pattern Recognition (CVPR), 2015
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_nakano_icpr_2025_one_sided",
        &lomsac_nakano_icpr_2025_one_sided_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] Gaku Nakano. "Inverse DLT Method for One-Sided Radial Distortion Homography", In
                International Conference on Pattern Recognition (ICPR), 2024.
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_wadenback_3dv_2026_one_sided",
        &lomsac_wadenback_3dv_2026_one_sided_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] TODO
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );

    m.def(
        "lomsac_wadenback_3dv_2026_two_sided_equal",
        &lomsac_wadenback_3dv_2026_two_sided_equal_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] TODO
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
    m.def(
        "lomsac_wadenback_3dv_2026_two_sided",
        &lomsac_wadenback_3dv_2026_two_sided_wrapper,
        R"pbdoc(
            Solver from [1] in a LOMSAC framework [2].
            
            [1] TODO
            [2] Karel Lebeda, Jiri Matas, and Ondrej Chum. "Fixing the Locally Optimized RANSAC", In the
            Proceedings of the British Machine Vision Conference (BMVC), 2012.
        )pbdoc",
        "x"_a,
        "y"_a,
        "options"_a
    );
};
