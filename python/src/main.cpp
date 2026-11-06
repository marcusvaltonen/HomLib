#include "binding_utils.h"

#include <RansacLib/ransac.h>
#include <lomsac_runner.h>
#include <generate_problem_instance.hpp>

using namespace pybind11::literals;

// --- Forward declarations for backend binding functions ---

#ifdef HOMLIB_HAS_SUPERANSAC
void bind_superansac(py::module_& m);
#endif
#ifdef HOMLIB_HAS_GCRANSAC
void bind_gcransac(py::module_& m);
#endif
#ifdef HOMLIB_HAS_MAGSAC
void bind_magsac(py::module_& m);
#endif

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
        auto s = make_solver<Solver>(solver);
        std::vector<HomLib::PoseData> poses;
        s.solve(x, y, &poses, direct_range(solver, x));
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
        auto s = make_solver<Solver>(solver);
        std::vector<HomLib::PoseData> poses;
        s.solve(x, y, A, &poses, direct_range(solver, x));
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
        auto s = make_solver<Solver>(solver);
        std::vector<HomLib::PoseData> poses;
        s.solve(x, y, ori, &poses, direct_range(solver, x));
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
        auto s = make_solver<Solver>(solver);
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
        auto s = make_solver<Solver>(solver);
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
        auto s = make_solver<Solver>(solver);
        return HomLib::run_lomsac_orientation(s, x, y, ori, options);
    });
}


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
        .value("CaiPAMI2025", HomLib::SolverId::CaiPAMI2025)
        .value("DLT", HomLib::SolverId::DLT)
        .value("FitzgibbonCVPR2001OneSided", HomLib::SolverId::FitzgibbonCVPR2001OneSided)
        .value("FitzgibbonCVPR2001TwoSidedEqual", HomLib::SolverId::FitzgibbonCVPR2001TwoSidedEqual)
        .value("KukelovaCVPR2015TwoSidedEqual", HomLib::SolverId::KukelovaCVPR2015TwoSidedEqual)
        .value("KukelovaCVPR2015TwoSidedEqual6Pt", HomLib::SolverId::KukelovaCVPR2015TwoSidedEqual6Pt)
        .value("KukelovaCVPR2015TwoSided", HomLib::SolverId::KukelovaCVPR2015TwoSided)
        .value("KukelovaCVPR2015TwoSided6Pt", HomLib::SolverId::KukelovaCVPR2015TwoSided6Pt)
        .value("NakanoICPR2025OneSided", HomLib::SolverId::NakanoICPR2025OneSided)
        .value("NakanoICPR2025OneSidedRight", HomLib::SolverId::NakanoICPR2025OneSidedRight)
        .value("FitzgibbonCVPR2001OneSidedRight", HomLib::SolverId::FitzgibbonCVPR2001OneSidedRight)
        .value("Wadenback3DV2026OneSided", HomLib::SolverId::Wadenback3DV2026OneSided)
        .value("Wadenback3DV2026OneSidedRight", HomLib::SolverId::Wadenback3DV2026OneSidedRight)
        .value("Wadenback3DV2026TwoSidedEqual", HomLib::SolverId::Wadenback3DV2026TwoSidedEqual)
        .value("Wadenback3DV2026TwoSided", HomLib::SolverId::Wadenback3DV2026TwoSided)
        .value("GuoPRL2019", HomLib::SolverId::GuoPRL2019)
        .value("BarathVISAPP2016Affine", HomLib::SolverId::BarathVISAPP2016Affine)
        .value("ValtonenOrnhagICPR2026AffineSingleSided", HomLib::SolverId::ValtonenOrnhagICPR2026AffineSingleSided)
        .value("ValtonenOrnhagICPR2026AffineSingleSidedLeft", HomLib::SolverId::ValtonenOrnhagICPR2026AffineSingleSidedLeft)
        .value("ValtonenOrnhagICPR2026OrientationSingleSided", HomLib::SolverId::ValtonenOrnhagICPR2026OrientationSingleSided)
        .value("ValtonenOrnhagICPR2026OrientationSingleSidedLeft", HomLib::SolverId::ValtonenOrnhagICPR2026OrientationSingleSidedLeft);

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
        .def(py::init<HomLib::SolverId, py::kwargs>(),
            R"pbdoc(
                Construct a Solver.

                Parameters
                ----------
                id : SolverId
                    Which solver algorithm to use.
                **kwargs
                    Solver-specific options (e.g. ``extra_check=True``).
                    Unknown options are silently ignored.
            )pbdoc",
            "id"_a)
        .def_readonly("id", &PySolver::id, "Solver identifier (:class:`SolverId`).")
        .def_readonly("feature_type", &PySolver::feature_type, "Feature type (:class:`FeatureType`).")
        .def_readonly("options", &PySolver::options, "Solver-specific options.")
        .def_readwrite("normalize_image_coord", &PySolver::normalize_image_coord,
            R"pbdoc(
                Whether the robust estimators normalize each sample (default True).

                They divide the points of both images by ``mean(|x|) / sqrt(2)``, computed from the
                first image's points (measured from the distortion center, the origin), and the
                range on the distortion (:meth:`set_k_range`) applies in these normalized
                coordinates. :meth:`estimate` doesn't normalize, but converts a range set with
                :meth:`set_k_range` to the points' scale the same way.

                With False, the points are used as given, and the range applies in their
                coordinates, in :meth:`estimate` too. The default [-1, 0.1] then rejects almost
                nothing for pixel coordinates, so set a range suited to them. This is useful when
                each image's focal length is known: with each image's centered points divided by
                it, a range chosen for the lens applies directly, and images with very different
                pixel scales (e.g. a large panorama and a small view of it) keep each image's k at
                its own scale. See the solver guide, "Range check on the distortion".
            )pbdoc")
        .def("set_k_range", [](PySolver& s, std::optional<double> k_min, std::optional<double> k_max) {
                const double lo = k_min ? *k_min : -std::numeric_limits<double>::infinity();
                const double hi = k_max ? *k_max : std::numeric_limits<double>::infinity();
                if (!(lo <= hi)) {
                    throw std::invalid_argument("k_min should be at most k_max (and neither NaN).");
                }
                s.has_k_range = true;
                s.k_min = lo;
                s.k_max = hi;
            },
            R"pbdoc(
                Keep only the minimal solutions with distortion parameters in [k_min, k_max].

                The bounds are in HomLib's normalized coordinates: the robust estimators (LOMSAC
                and the other backends) divide each sample by ``mean(|x|) / sqrt(2)``, computed
                from the first image's points, whatever the units of the input. They apply
                [-1, 0.1] by default; this replaces it. :meth:`estimate` calls the solver on the
                points as given, with no range unless one is set, converted to the points' scale
                (see :attr:`normalize_image_coord`). A parameter equal to 0, the one that a
                one-sided solver doesn't estimate, passes. The check selects the hypotheses of the
                minimal solver; it doesn't constrain the refinement.

                The default [-1, 0.1] only rejects roots that no real lens gives (see the solver
                guide, "Range check on the distortion"). ``set_k_range()`` without bounds, like
                infinite bounds on both sides, disables the check.

                Parameters
                ----------
                k_min : float or None
                    Lower bound (None: no lower bound).
                k_max : float or None
                    Upper bound (None: no upper bound). ``set_k_range()`` with no bounds disables
                    the check.
            )pbdoc",
            "k_min"_a = std::nullopt, "k_max"_a = std::nullopt)
        .def("reset_k_range", [](PySolver& s) {
                s.has_k_range = false;
                s.k_min = -std::numeric_limits<double>::infinity();
                s.k_max = std::numeric_limits<double>::infinity();
            },
            R"pbdoc(
                Back to the default: [-1, 0.1] in the robust estimators, no range in
                :meth:`estimate` (see :meth:`set_k_range`).
            )pbdoc")
        .def_property_readonly("k_range", [](const PySolver& s) -> py::object {
                if (!s.has_k_range) {
                    return py::none();
                }
                return py::make_tuple(s.k_min, s.k_max);
            },
            "The range set by :meth:`set_k_range`, as (k_min, k_max), or None if none was set (then [-1, 0.1] in the robust "
            "estimators, no range in :meth:`estimate`).")
        .def("__repr__", [](const PySolver& s) {
            auto py_id = py::cast(s.id);
            std::string name = py_id.attr("name").cast<std::string>();
            std::string r = "Solver(SolverId." + name;
            for (auto item : s.options) {
                r += ", " + py::str(item.first).cast<std::string>() +
                     "=" + py::str(item.second).cast<std::string>();
            }
            if (s.has_k_range) {
                r += ", k_range=" + py::repr(py::make_tuple(s.k_min, s.k_max)).cast<std::string>();
            }
            if (!s.normalize_image_coord) {
                r += ", normalize_image_coord=False";
            }
            return r + ")";
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
        .def_readwrite("ori", &HomLib::ProblemInstance::ori, "Orientation vectors (2D).")
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
        .def_property_readonly("inliers",
            [](const ransac_lib::RansacStatistics &s) {
                return static_cast<int>(s.inlier_indices.size());
            },
            "Number of inliers (convenience for len(inlier_indices)).")
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

    // --- Backend bindings (each in its own compilation unit) ---

#ifdef HOMLIB_HAS_SUPERANSAC
    bind_superansac(m);
#else
    m.attr("_has_superansac") = false;
#endif

#ifdef HOMLIB_HAS_GCRANSAC
    bind_gcransac(m);
#else
    m.attr("_has_gcransac") = false;
#endif

#ifdef HOMLIB_HAS_MAGSAC
    bind_magsac(m);
#else
    m.attr("_has_magsac") = false;
#endif
};
