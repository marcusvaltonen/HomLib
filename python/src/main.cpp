#include <tuple>
#include <vector>
#include <sstream>
#include <type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include <Eigen/Dense>

#include <RansacLib/ransac.h>

#include <posedata.hpp>
#include <solver_id.h>
#include <solver_dispatch.h>
#include <lomsac_runner.h>
#include <generate_problem_instance.hpp>

#ifdef HOMLIB_HAS_SUPERANSAC
#include <superansac_runner.h>
#endif

namespace py = pybind11;
using namespace pybind11::literals;

// --- Helpers ---

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

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

void preprocess_affine(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic> &A_,
    std::vector<Eigen::Vector2d> *x,
    std::vector<Eigen::Vector2d> *y,
    std::vector<Eigen::Matrix2d> *A
) {
    if (x_.cols() != y_.cols()) {
        throw std::invalid_argument("x and y should be of equal size.");
    }
    if (x_.cols() != A_.cols()) {
        throw std::invalid_argument("x and A should be of equal size.");
    }
    for (size_t i=0; i < x_.cols(); i++) {
        x->push_back(x_.col(i));
        y->push_back(y_.col(i));
        Eigen::Matrix2d tmp;
        tmp << A_.col(i);
        tmp.transposeInPlace();
        A->push_back(tmp);
    }
}

void preprocess_orientation(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &ori_,
    std::vector<Eigen::Vector2d> *x,
    std::vector<Eigen::Vector2d> *y,
    std::vector<Eigen::Vector2d> *ori
) {
    if (x_.cols() != y_.cols()) {
        throw std::invalid_argument("x and y should be of equal size.");
    }
    if (x_.cols() != ori_.cols()) {
        throw std::invalid_argument("x and ori should be of equal size.");
    }
    for (size_t i=0; i < x_.cols(); i++) {
        x->push_back(x_.col(i));
        y->push_back(y_.col(i));
        ori->push_back(ori_.col(i));
    }
}

// --- Solver construction with optional extra_check ---

template <typename Solver>
Solver make_solver(bool extra_check) {
    if constexpr (std::is_constructible_v<Solver, bool>) {
        return Solver(extra_check);
    } else {
        return Solver();
    }
}

// --- PySolver: wraps SolverId + params for use from Python ---

struct PySolver {
    HomLib::SolverId id;
    HomLib::FeatureType feature_type;
    bool extra_check;

    PySolver(HomLib::SolverId id, bool extra_check = false)
        : id(id), feature_type(HomLib::feature_type_for(id)), extra_check(extra_check) {}
};

// --- Dispatch: direct estimation ---

std::vector<HomLib::PoseData> dispatch_estimate_point(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_)
{
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);

    return HomLib::dispatch_point_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        std::vector<HomLib::PoseData> poses;
        s.solve(x, y, &poses);
        return poses;
    });
}

std::vector<HomLib::PoseData> dispatch_estimate_affine(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic>& A_)
{
    std::vector<Eigen::Vector2d> x, y;
    std::vector<Eigen::Matrix2d> A;
    preprocess_affine(x_, y_, A_, &x, &y, &A);

    return HomLib::dispatch_affine_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        std::vector<HomLib::PoseData> poses;
        s.solve(x, y, A, &poses);
        return poses;
    });
}

std::vector<HomLib::PoseData> dispatch_estimate_orientation(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& ori_)
{
    std::vector<Eigen::Vector2d> x, y, ori;
    preprocess_orientation(x_, y_, ori_, &x, &y, &ori);

    return HomLib::dispatch_orientation_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        std::vector<HomLib::PoseData> poses;
        s.solve(x, y, ori, &poses);
        return poses;
    });
}

// --- Dispatch: LOMSAC ---

std::pair<HomLib::PoseData, ransac_lib::RansacStatistics> dispatch_lomsac_point(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const ransac_lib::LORansacOptions& options)
{
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);

    return HomLib::dispatch_point_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_lomsac_point(s, x, y, options);
    });
}

std::pair<HomLib::PoseData, ransac_lib::RansacStatistics> dispatch_lomsac_affine(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic>& A_,
    const ransac_lib::LORansacOptions& options)
{
    std::vector<Eigen::Vector2d> x, y;
    std::vector<Eigen::Matrix2d> A;
    preprocess_affine(x_, y_, A_, &x, &y, &A);

    return HomLib::dispatch_affine_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_lomsac_affine(s, x, y, A, options);
    });
}

std::pair<HomLib::PoseData, ransac_lib::RansacStatistics> dispatch_lomsac_orientation(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& ori_,
    const ransac_lib::LORansacOptions& options)
{
    std::vector<Eigen::Vector2d> x, y, ori;
    preprocess_orientation(x_, y_, ori_, &x, &y, &ori);

    return HomLib::dispatch_orientation_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_lomsac_orientation(s, x, y, ori, options);
    });
}


// --- Dispatch: SupeRANSAC (conditional) ---

#ifdef HOMLIB_HAS_SUPERANSAC
std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_point(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    preprocess(x_, y_, &x, &y);

    return HomLib::dispatch_point_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_superansac_point(s, x, y, settings);
    });
}

std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_affine(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 4, Eigen::Dynamic>& A_,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y;
    std::vector<Eigen::Matrix2d> A;
    preprocess_affine(x_, y_, A_, &x, &y, &A);

    return HomLib::dispatch_affine_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_superansac_affine(s, x, y, A, settings);
    });
}

std::pair<HomLib::PoseData, HomLib::SuperansacStatistics> dispatch_superansac_orientation(
    const PySolver& solver,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& x_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& y_,
    const Eigen::Matrix<double, 2, Eigen::Dynamic>& ori_,
    const superansac::RANSACSettings& settings)
{
    std::vector<Eigen::Vector2d> x, y, ori;
    preprocess_orientation(x_, y_, ori_, &x, &y, &ori);

    return HomLib::dispatch_orientation_solver(solver.id, [&](auto tag) {
        using Solver = std::decay_t<decltype(tag)>;
        auto s = make_solver<Solver>(solver.extra_check);
        return HomLib::run_superansac_orientation(s, x, y, ori, settings);
    });
}
#endif // HOMLIB_HAS_SUPERANSAC


PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Wrappers for C++ core algorithms of HomLib.
    )pbdoc";

    #ifdef HOMLIB_VERSION
        m.attr("__version__") = py::str(HOMLIB_VERSION);
    #else
        m.attr("__version__") = py::str("dev");
    #endif

    // --- Enums ---

    py::enum_<HomLib::SolverId>(m, "SolverId",
        R"pbdoc(
            Identifier for a specific minimal solver algorithm.
        )pbdoc")
        .value("FitzgibbonCVPR2001OneSided", HomLib::SolverId::FitzgibbonCVPR2001OneSided)
        .value("FitzgibbonCVPR2001TwoSidedEqual", HomLib::SolverId::FitzgibbonCVPR2001TwoSidedEqual)
        .value("KukelovaCVPR2015TwoSidedEqual", HomLib::SolverId::KukelovaCVPR2015TwoSidedEqual)
        .value("KukelovaCVPR2015TwoSidedEqual6Pt", HomLib::SolverId::KukelovaCVPR2015TwoSidedEqual6Pt)
        .value("KukelovaCVPR2015TwoSided", HomLib::SolverId::KukelovaCVPR2015TwoSided)
        .value("KukelovaCVPR2015TwoSided6Pt", HomLib::SolverId::KukelovaCVPR2015TwoSided6Pt)
        .value("NakanoICPR2025OneSided", HomLib::SolverId::NakanoICPR2025OneSided)
        .value("NakanoICPR2025OneSidedRight", HomLib::SolverId::NakanoICPR2025OneSidedRight)
        .value("Wadenback3DV2026OneSided", HomLib::SolverId::Wadenback3DV2026OneSided)
        .value("Wadenback3DV2026TwoSidedEqual", HomLib::SolverId::Wadenback3DV2026TwoSidedEqual)
        .value("Wadenback3DV2026TwoSided", HomLib::SolverId::Wadenback3DV2026TwoSided)
        .value("BarathVISAPP2016Affine", HomLib::SolverId::BarathVISAPP2016Affine)
        .value("ValtonenOrnhagICPR2026AffineSingleSided", HomLib::SolverId::ValtonenOrnhagICPR2026AffineSingleSided)
        .value("ValtonenOrnhagICPR2026OrientationSingleSided", HomLib::SolverId::ValtonenOrnhagICPR2026OrientationSingleSided);

    py::enum_<HomLib::FeatureType>(m, "FeatureType",
        R"pbdoc(
            Type of input features expected by a solver.
        )pbdoc")
        .value("Points", HomLib::FeatureType::Points)
        .value("Affine", HomLib::FeatureType::Affine)
        .value("Orientation", HomLib::FeatureType::Orientation);

    // --- PySolver ---

    py::class_<PySolver>(m, "Solver",
        R"pbdoc(
            A solver for homography estimation with optional radial distortion correction.

            Wraps a specific solver algorithm identified by :class:`SolverId`.
            Solver-specific parameters (e.g. ``extra_check``) are set at construction.
        )pbdoc")
        .def(py::init<HomLib::SolverId, bool>(),
            R"pbdoc(
                Construct a Solver.

                Parameters
                ----------
                id : SolverId
                    Which solver algorithm to use.
                extra_check : bool
                    Whether to apply extra validation checks (only affects solvers that support it).
            )pbdoc",
            "id"_a, "extra_check"_a = false)
        .def_readonly("id", &PySolver::id, "Solver identifier (:class:`SolverId`).")
        .def_readonly("feature_type", &PySolver::feature_type, "Feature type (:class:`FeatureType`).")
        .def_readonly("extra_check", &PySolver::extra_check, "Whether extra validation checks are enabled.")
        .def("__repr__", [](const PySolver& s) {
            auto py_id = py::cast(s.id);
            std::string name = py_id.attr("name").cast<std::string>();
            return "Solver(SolverId." + name +
                   ", extra_check=" + std::string(s.extra_check ? "True" : "False") + ")";
        });

    // --- PoseData ---

    py::class_<HomLib::PoseData>(m, "PoseData",
        R"pbdoc(
            Solver output containing the estimated homography and distortion parameters.
        )pbdoc")
        .def(
            py::init<Eigen::Matrix3d, double, double, double>(),
            R"pbdoc(
                Construct a PoseData instance.

                Parameters
                ----------
                homography : numpy.ndarray
                    3x3 homography matrix.
                focal_length : float
                    Focal length.
                distortion_parameter : float
                    Distortion parameter for the first image.
                distortion_parameter2 : float
                    Distortion parameter for the second image.
            )pbdoc",
            "homography"_a,
            "focal_length"_a,
            "distortion_parameter"_a,
            "distortion_parameter2"_a
        )
        .def_readwrite("homography", &HomLib::PoseData::homography, "3x3 homography matrix (numpy.ndarray).")
        .def_readwrite("focal_length", &HomLib::PoseData::focal_length, "Estimated focal length.")
        .def_readwrite("distortion_parameter", &HomLib::PoseData::distortion_parameter,
            "Radial distortion parameter for the first set of points (x).")
        .def_readwrite("distortion_parameter2", &HomLib::PoseData::distortion_parameter2,
            "Radial distortion parameter for the second set of points (y).")
        .def("__repr__",
            [](const HomLib::PoseData &p) {
                return "PoseData("
                    "H=[3x3 np.array], "
                    "focal_length=" + std::to_string(p.focal_length) + ", "
                    "distortion_parameter=" + to_string_with_precision(p.distortion_parameter, 16) + ", "
                    "distortion_parameter2=" + to_string_with_precision(p.distortion_parameter2, 16) +
                    ")";
            }
        );

    // --- DistortionCase ---

    py::enum_<HomLib::DistortionCase>(m, "DistortionCase",
        R"pbdoc(
            Radial distortion model used when generating synthetic problems.

            "Left" and "right" refer to the side of the homography equation
            ``undistort(y) = H * x``, not to image ordering. Left-sided means
            distortion is on y (x2); right-sided means distortion is on x (x1).
        )pbdoc")
        .value("NO_DISTORTION", HomLib::DistortionCase::NO_DISTORTION, "No radial distortion.")
        .value("ONE_SIDED_LEFT", HomLib::DistortionCase::ONE_SIDED_LEFT,
            "Distortion on x2 (left-hand side of the equation); x1 is undistorted.")
        .value("ONE_SIDED_RIGHT", HomLib::DistortionCase::ONE_SIDED_RIGHT,
            "Distortion on x1 (right-hand side of the equation); x2 is undistorted.")
        .value("TWO_SIDED_EQUAL", HomLib::DistortionCase::TWO_SIDED_EQUAL,
            "Equal distortion on both sets of points.")
        .value("TWO_SIDED", HomLib::DistortionCase::TWO_SIDED,
            "Independent distortion on each set of points.")
        .export_values();

    // --- ProblemConfig ---

    py::class_<HomLib::ProblemConfig>(m, "ProblemConfig",
        R"pbdoc(
            Configuration for synthetic problem generation.
        )pbdoc")
        .def(
            py::init<HomLib::DistortionCase, double, int>(),
            R"pbdoc(
                Construct a ProblemConfig.

                Parameters
                ----------
                distortion : DistortionCase
                    Which distortion model to use.
                point_noise : float
                    Standard deviation of Gaussian noise on point coordinates.
                number_points : int
                    Number of point correspondences.
            )pbdoc",
            "distortion"_a,
            "point_noise"_a,
            "number_points"_a
        )
        .def_readwrite("distortion", &HomLib::ProblemConfig::distortion, "Distortion model (:class:`DistortionCase`).")
        .def_readwrite("point_noise", &HomLib::ProblemConfig::point_noise, "Gaussian noise standard deviation.")
        .def_readwrite("number_points", &HomLib::ProblemConfig::number_points, "Number of point correspondences.")
        .def_readwrite("camera_fov", &HomLib::ProblemConfig::camera_fov_, "Camera field of view in degrees.")
        .def_readwrite("min_depth", &HomLib::ProblemConfig::min_depth_, "Minimum scene depth.")
        .def_readwrite("max_depth", &HomLib::ProblemConfig::max_depth_, "Maximum scene depth.")
        .def_readwrite("min_focal", &HomLib::ProblemConfig::min_focal_, "Minimum focal length.")
        .def_readwrite("max_focal", &HomLib::ProblemConfig::max_focal_, "Maximum focal length.")
        .def_readwrite("min_dist", &HomLib::ProblemConfig::min_dist_, "Minimum distortion parameter.")
        .def_readwrite("max_dist", &HomLib::ProblemConfig::max_dist_, "Maximum distortion parameter.")
        .def("__repr__",
            [](const HomLib::ProblemConfig &p) {
                std::string distortion;
                switch (p.distortion) {
            	    case HomLib::DistortionCase::NO_DISTORTION:
                        distortion = "NO_DISTORTION";
		                break;
            	    case HomLib::DistortionCase::ONE_SIDED_LEFT:
                        distortion = "ONE_SIDED_LEFT";
		                break;
            	    case HomLib::DistortionCase::ONE_SIDED_RIGHT:
                        distortion = "ONE_SIDED_RIGHT";
		                break;
            	    case HomLib::DistortionCase::TWO_SIDED_EQUAL:
                        distortion = "TWO_SIDED_EQUAL";
		                break;
            	    case HomLib::DistortionCase::TWO_SIDED:
            	        distortion = "TWO_SIDED";
		                break;
        	    }
                return "ProblemConfig("
                    "distortion=" + distortion + ", "
                    "point_noise=" + to_string_with_precision(p.point_noise) + ", "
                    "number_points=" +std::to_string(p.number_points) + ", "
                    "camera_fov=" + to_string_with_precision(p.camera_fov_, 2) + ", "
                    "min_depth=" + to_string_with_precision(p.min_depth_, 2) + ", "
                    "max_depth=" + to_string_with_precision(p.max_depth_, 2) + ", "
                    "min_focal=" + to_string_with_precision(p.min_focal_, 2) + ", "
                    "max_focal=" + to_string_with_precision(p.max_focal_, 2) + ", "
                    "min_dist=" + to_string_with_precision(p.min_dist_, 2) + ", "
                    "max_dist=" + to_string_with_precision(p.max_dist_, 2) +
                    ")";
            }
        );

    // --- ProblemInstance ---

    py::class_<HomLib::ProblemInstance>(m, "ProblemInstance",
        R"pbdoc(
            A synthetic problem instance with ground-truth data and point correspondences.
        )pbdoc")
        .def(
            py::init<HomLib::PoseData, std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector2d>, std::vector<Eigen::Matrix2d>>(),
            "Construct a ProblemInstance.",
            "posedata"_a,
            "x1"_a,
            "x2"_a,
            "A"_a
        )
        .def_readwrite("posedata", &HomLib::ProblemInstance::posedata, "Ground-truth :class:`PoseData`.")
        .def_readwrite("x1", &HomLib::ProblemInstance::x1, "First set of 2D points.")
        .def_readwrite("x2", &HomLib::ProblemInstance::x2, "Second set of 2D points.")
        .def_readwrite("A", &HomLib::ProblemInstance::A, "Local affine transformations (2x2).")
        .def("hom_error", &HomLib::ProblemInstance::hom_error,
            R"pbdoc(
                Homography error between an estimate and the ground truth.

                Parameters
                ----------
                H_est : numpy.ndarray
                    Estimated 3x3 homography matrix.

                Returns
                -------
                float
                    Frobenius norm between the normalised ground-truth and estimated homographies.
            )pbdoc",
            "H_est"_a)
        .def("dist_error", &HomLib::ProblemInstance::dist_error,
            R"pbdoc(
                Mean absolute distortion error.

                Parameters
                ----------
                k1_est : float
                    Estimated distortion parameter for the first image.
                k2_est : float
                    Estimated distortion parameter for the second image.

                Returns
                -------
                float
                    ``0.5 * (|k1_gt - k1_est| + |k2_gt - k2_est|)``.
            )pbdoc",
            "k1_est"_a, "k2_est"_a)
        .def("__repr__",
            [](const HomLib::ProblemInstance &p) {
                return "ProblemInstance("
                    "posedata=PoseData(), "
                    "x1=[list of length " + std::to_string(p.x1.size()) + "], "
                    "x2=[list of length " + std::to_string(p.x2.size()) + "], "
                    "A=[list of length " + std::to_string(p.A.size()) + "]"
                    ")";
            }
        );
    m.def(
        "generate_problem_instance",
        &HomLib::generate_problem_instance,
        R"pbdoc(
            Generate a synthetic problem instance with a random homography and point correspondences.

            Creates a random 3D plane, projects points through two cameras, and applies
            radial distortion according to the configuration. The returned instance
            contains ground-truth data and can be used for solver evaluation.

            Parameters
            ----------
            settings : ProblemConfig
                Configuration controlling distortion model, noise level, and number of points.

            Returns
            -------
            ProblemInstance
                The generated problem instance with ground-truth ``posedata``, distorted
                point correspondences ``x1`` and ``x2``, and local affine data ``A``.
        )pbdoc",
        "settings"_a
    );

    // --- LORansacOptions ---

    py::class_<ransac_lib::LORansacOptions>(m, "LORansacOptions",
        R"pbdoc(
            Options for the Locally Optimized MSAC (LOMSAC) estimator.
        )pbdoc")
        .def(
            py::init<>()
        )
        .def_readwrite("min_num_iterations", &ransac_lib::LORansacOptions::min_num_iterations_,
            "Minimum number of RANSAC iterations.")
        .def_readwrite("max_num_iterations", &ransac_lib::LORansacOptions::max_num_iterations_,
            "Maximum number of RANSAC iterations.")
        .def_readwrite("success_probability", &ransac_lib::LORansacOptions::success_probability_,
            "Desired probability of finding the correct model.")
        .def_readwrite("squared_inlier_threshold", &ransac_lib::LORansacOptions::squared_inlier_threshold_,
            "Squared inlier threshold for the MSAC score.")
        .def_readwrite("random_seed", &ransac_lib::LORansacOptions::random_seed_,
            "Seed for the random number generator.")
        .def_readwrite("num_lo_steps", &ransac_lib::LORansacOptions::num_lo_steps_,
            "Number of local optimisation steps.")
        .def_readwrite("threshold_multiplier", &ransac_lib::LORansacOptions::threshold_multiplier_,
            "Multiplier applied to the threshold during local optimisation.")
        .def_readwrite("num_lsq_iterations", &ransac_lib::LORansacOptions::num_lsq_iterations_,
            "Number of least-squares iterations in local optimisation.")
        .def_readwrite("min_sample_multiplicator", &ransac_lib::LORansacOptions::min_sample_multiplicator_,
            "Multiplicator for the minimum sample size.")
        .def_readwrite("non_min_sample_multiplier", &ransac_lib::LORansacOptions::non_min_sample_multiplier_,
            "Multiplier for the non-minimal sample size in local optimisation.")
        .def_readwrite("lo_starting_iterations", &ransac_lib::LORansacOptions::lo_starting_iterations_,
            "Number of initial iterations before local optimisation begins.")
        .def_readwrite("final_least_squares", &ransac_lib::LORansacOptions::final_least_squares_,
            "Whether to run a final least-squares refinement.")
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

    // --- RansacStatistics ---

    py::class_<ransac_lib::RansacStatistics>(m, "RansacStatistics",
        R"pbdoc(
            Statistics from a RANSAC / LOMSAC estimation run.
        )pbdoc")
        .def(
            py::init<>()
        )
        .def_readwrite("num_iterations", &ransac_lib::RansacStatistics::num_iterations,
            "Total number of RANSAC iterations.")
        .def_readwrite("best_num_inliers", &ransac_lib::RansacStatistics::best_num_inliers,
            "Number of inliers for the best model.")
        .def_readwrite("best_model_score", &ransac_lib::RansacStatistics::best_model_score,
            "MSAC score of the best model.")
        .def_readwrite("inlier_ratio", &ransac_lib::RansacStatistics::inlier_ratio,
            "Ratio of inliers to total correspondences.")
        .def_readwrite("inlier_indices", &ransac_lib::RansacStatistics::inlier_indices,
            "Indices of inlier correspondences.")
        .def_readwrite("number_lo_iterations", &ransac_lib::RansacStatistics::number_lo_iterations,
            "Number of local optimisation iterations.")
        .def("__repr__",
            [](const ransac_lib::RansacStatistics &s) {
                return "RansacStatistics("
                    "num_iterations=" + std::to_string(s.num_iterations) + ", "
                    "best_num_inliers=" + std::to_string(s.best_num_inliers) + ", "
                    "best_model_score=" + std::to_string(s.best_model_score) + ", "
                    "inlier_ratio=" + std::to_string(s.inlier_ratio) + ", "
                    "inlier_indices=[int list of length " + std::to_string(s.inlier_indices.size()) + "], "
                    "number_lo_iterations=" + std::to_string(s.number_lo_iterations) + ", "
                    ")";
            }
        );

    // --- Dispatch functions ---

    m.def("dispatch_estimate_point", &dispatch_estimate_point,
        R"pbdoc(
            Run a point-based minimal solver directly (no RANSAC).

            Parameters
            ----------
            solver : Solver
                A point-based solver instance.
            x : numpy.ndarray
                2xN array of points in the first image.
            y : numpy.ndarray
                2xN array of points in the second image.

            Returns
            -------
            list[PoseData]
                List of candidate solutions.
        )pbdoc",
        "solver"_a, "x"_a, "y"_a);

    m.def("dispatch_estimate_affine", &dispatch_estimate_affine,
        R"pbdoc(
            Run an affine-covariant minimal solver directly (no RANSAC).

            Parameters
            ----------
            solver : Solver
                An affine-covariant solver instance.
            x : numpy.ndarray
                2xN array of points in the first image.
            y : numpy.ndarray
                2xN array of points in the second image.
            A : numpy.ndarray
                4xN array of local affine transformations (row-major 2x2 per column).

            Returns
            -------
            list[PoseData]
                List of candidate solutions.
        )pbdoc",
        "solver"_a, "x"_a, "y"_a, "A"_a);

    m.def("dispatch_estimate_orientation", &dispatch_estimate_orientation,
        R"pbdoc(
            Run an orientation-covariant minimal solver directly (no RANSAC).

            Parameters
            ----------
            solver : Solver
                An orientation-covariant solver instance.
            x : numpy.ndarray
                2xN array of points in the first image.
            y : numpy.ndarray
                2xN array of points in the second image.
            ori : numpy.ndarray
                2xN array of orientation vectors.

            Returns
            -------
            list[PoseData]
                List of candidate solutions.
        )pbdoc",
        "solver"_a, "x"_a, "y"_a, "ori"_a);

    m.def("dispatch_lomsac_point", &dispatch_lomsac_point,
        R"pbdoc(
            Run a point-based solver within the LOMSAC framework.

            Parameters
            ----------
            solver : Solver
                A point-based solver instance.
            x : numpy.ndarray
                2xN array of points in the first image.
            y : numpy.ndarray
                2xN array of points in the second image.
            options : LORansacOptions
                LOMSAC options.

            Returns
            -------
            tuple[PoseData, RansacStatistics]
                Best model and estimation statistics.
        )pbdoc",
        "solver"_a, "x"_a, "y"_a, "options"_a);

    m.def("dispatch_lomsac_affine", &dispatch_lomsac_affine,
        R"pbdoc(
            Run an affine-covariant solver within the LOMSAC framework.

            Parameters
            ----------
            solver : Solver
                An affine-covariant solver instance.
            x : numpy.ndarray
                2xN array of points in the first image.
            y : numpy.ndarray
                2xN array of points in the second image.
            A : numpy.ndarray
                4xN array of local affine transformations.
            options : LORansacOptions
                LOMSAC options.

            Returns
            -------
            tuple[PoseData, RansacStatistics]
                Best model and estimation statistics.
        )pbdoc",
        "solver"_a, "x"_a, "y"_a, "A"_a, "options"_a);

    m.def("dispatch_lomsac_orientation", &dispatch_lomsac_orientation,
        R"pbdoc(
            Run an orientation-covariant solver within the LOMSAC framework.

            Parameters
            ----------
            solver : Solver
                An orientation-covariant solver instance.
            x : numpy.ndarray
                2xN array of points in the first image.
            y : numpy.ndarray
                2xN array of points in the second image.
            ori : numpy.ndarray
                2xN array of orientation vectors.
            options : LORansacOptions
                LOMSAC options.

            Returns
            -------
            tuple[PoseData, RansacStatistics]
                Best model and estimation statistics.
        )pbdoc",
        "solver"_a, "x"_a, "y"_a, "ori"_a, "options"_a);

    // --- SupeRANSAC bindings (conditional) ---

#ifdef HOMLIB_HAS_SUPERANSAC
    m.attr("_has_superansac") = true;

    py::class_<HomLib::SuperansacStatistics>(m, "SuperansacStatistics",
        R"pbdoc(
            Statistics from a SupeRANSAC estimation run.
        )pbdoc")
        .def(py::init<>())
        .def_readwrite("num_iterations", &HomLib::SuperansacStatistics::num_iterations,
            "Total number of RANSAC iterations.")
        .def_readwrite("inliers", &HomLib::SuperansacStatistics::inliers,
            "Indices of inlier correspondences.")
        .def_readwrite("best_score", &HomLib::SuperansacStatistics::best_score,
            "Score of the best model.")
        .def("__repr__",
            [](const HomLib::SuperansacStatistics& s) {
                return "SuperansacStatistics("
                    "num_iterations=" + std::to_string(s.num_iterations) + ", "
                    "inliers=[list of length " + std::to_string(s.inliers.size()) + "], "
                    "best_score=" + std::to_string(s.best_score) +
                    ")";
            });

    py::class_<superansac::RANSACSettings>(m, "SuperansacSettings",
        R"pbdoc(
            Settings for SupeRANSAC estimation.
        )pbdoc")
        .def(py::init<>())
        .def_readwrite("inlier_threshold", &superansac::RANSACSettings::inlierThreshold,
            "Inlier threshold.")
        .def_readwrite("confidence", &superansac::RANSACSettings::confidence,
            "Desired confidence level.")
        .def_readwrite("min_iterations", &superansac::RANSACSettings::minIterations,
            "Minimum number of iterations.")
        .def_readwrite("max_iterations", &superansac::RANSACSettings::maxIterations,
            "Maximum number of iterations.")
        .def_readwrite("top_k", &superansac::RANSACSettings::topKForLocalOptimization,
            "Number of top models for local optimization.")
        .def_readwrite("use_sprt", &superansac::RANSACSettings::useSprt,
            "Whether to use SPRT test.")
        .def_readwrite("sampler", &superansac::RANSACSettings::sampler,
            "Sampler type (default: PROSAC). Use SuperansacSamplerType enum.")
        .def_readwrite("scoring", &superansac::RANSACSettings::scoring,
            "Scoring type (default: MAGSAC). Use SuperansacScoringType enum.")
        .def_readwrite("local_optimization", &superansac::RANSACSettings::localOptimization,
            "Local optimization type. Use SuperansacLocalOptType enum.")
        .def_readwrite("final_optimization", &superansac::RANSACSettings::finalOptimization,
            "Final optimization type. Use SuperansacLocalOptType enum.")
        .def("__repr__",
            [](const superansac::RANSACSettings& s) {
                return "SuperansacSettings("
                    "inlier_threshold=" + std::to_string(s.inlierThreshold) + ", "
                    "confidence=" + std::to_string(s.confidence) + ", "
                    "min_iterations=" + std::to_string(s.minIterations) + ", "
                    "max_iterations=" + std::to_string(s.maxIterations) +
                    ")";
            });

    // Scoring type enum
    py::enum_<superansac::scoring::ScoringType>(m, "SuperansacScoringType")
        .value("RANSAC", superansac::scoring::ScoringType::RANSAC)
        .value("MSAC", superansac::scoring::ScoringType::MSAC)
        .value("MAGSAC", superansac::scoring::ScoringType::MAGSAC)
        .value("MINPRAN", superansac::scoring::ScoringType::MINPRAN)
        .value("ACRANSAC", superansac::scoring::ScoringType::ACRANSAC)
        .value("GAU", superansac::scoring::ScoringType::GAU)
        .value("ML", superansac::scoring::ScoringType::ML)
        .value("GRID", superansac::scoring::ScoringType::GRID);

    // Sampler type enum
    py::enum_<superansac::samplers::SamplerType>(m, "SuperansacSamplerType")
        .value("Uniform", superansac::samplers::SamplerType::Uniform)
        .value("PROSAC", superansac::samplers::SamplerType::PROSAC)
        .value("NAPSAC", superansac::samplers::SamplerType::NAPSAC)
        .value("ProgressiveNAPSAC", superansac::samplers::SamplerType::ProgressiveNAPSAC)
        .value("ImportanceSampler", superansac::samplers::SamplerType::ImportanceSampler)
        .value("ARSampler", superansac::samplers::SamplerType::ARSampler);

    // Local optimization type enum
    py::enum_<superansac::local_optimization::LocalOptimizationType>(m, "SuperansacLocalOptType")
        .value("None_", superansac::local_optimization::LocalOptimizationType::None)
        .value("LSQ", superansac::local_optimization::LocalOptimizationType::LSQ)
        .value("IRLS", superansac::local_optimization::LocalOptimizationType::IRLS)
        .value("NestedRANSAC", superansac::local_optimization::LocalOptimizationType::NestedRANSAC)
        .value("GCRANSAC", superansac::local_optimization::LocalOptimizationType::GCRANSAC)
        .value("IteratedLMEDS", superansac::local_optimization::LocalOptimizationType::IteratedLMEDS)
        .value("CrossValidation", superansac::local_optimization::LocalOptimizationType::CrossValidation);

    // Expose scoring/sampler/lo setters on settings
    py::class_<superansac::RANSACSettings> settings_class = m.attr("SuperansacSettings").cast<py::class_<superansac::RANSACSettings>>();

    m.def("dispatch_superansac_point", &dispatch_superansac_point,
        "Run a point-based solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "settings"_a);

    m.def("dispatch_superansac_affine", &dispatch_superansac_affine,
        "Run an affine-covariant solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "A"_a, "settings"_a);

    m.def("dispatch_superansac_orientation", &dispatch_superansac_orientation,
        "Run an orientation-covariant solver within the SupeRANSAC framework.",
        "solver"_a, "x"_a, "y"_a, "ori"_a, "settings"_a);
#else
    m.attr("_has_superansac") = false;
#endif // HOMLIB_HAS_SUPERANSAC
};
