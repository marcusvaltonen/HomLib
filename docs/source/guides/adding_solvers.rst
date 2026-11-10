.. _adding_solvers:

Adding a New Solver
===================

This guide covers the steps to add a new solver to HomLib, both for minimal
(RANSAC sample) and non-minimal (overdetermined) use.


Distortion-case base classes
----------------------------

Every point-based solver inherits from one of these intermediate base classes
in ``HomLib/robust/pose_estimator.h``:

.. list-table::
   :header-rows: 1
   :widths: 30 25 25

   * - Base class
     - ``DISTORTION_CASE``
     - Refinement
   * - ``NoDistortionPointEstimator<T>``
     - ``NO_DISTORTION``
     - ``refinement_no_dist``
   * - ``OneSidedLeftPointEstimator<T>``
     - ``ONE_SIDED_LEFT``
     - ``refinement_onesided``
   * - ``TwoSidedEqualPointEstimator<T>``
     - ``TWO_SIDED_EQUAL``
     - ``refinement_twosided_equal``
   * - ``TwoSidedPointEstimator<T>``
     - ``TWO_SIDED``
     - ``refinement_twosided``

Each base class provides:

- ``static constexpr DistortionCase DISTORTION_CASE`` — exposed to Python as
  ``Solver.distortion_case``.
- ``minimal_sample_size()`` — reads ``Solver::minimal_sample_size_value``.
- ``undistort()`` / ``distort()`` — identity or radial, matching the case.
- ``refine()`` — calls the correct Levenberg–Marquardt refinement function.

The solver only needs to implement ``solve()``.


Adding a minimal solver
-----------------------

1. **Create the solver header** ``HomLib/includes/HomLib/get_<name>.hpp``.
   Inherit from the appropriate base class:

   .. code-block:: cpp

      #include "pose_estimator.h"

      namespace HomLib {
      namespace MyPaper2026 {

      class SolverOneSided : public OneSidedLeftPointEstimator<SolverOneSided> {
      public:
          static constexpr int minimal_sample_size_value = 5;
          int solve(const std::vector<Eigen::Vector2d>& x,
                    const std::vector<Eigen::Vector2d>& y,
                    std::vector<HomLib::PoseData>* poses,
                    const HomLib::DistortionRange& range = {}) const {
              // Call your algorithm, push results to *poses
          }
      };

      } // namespace MyPaper2026
      } // namespace HomLib

2. **Add the enum value** in ``HomLib/robust/solver_id.h``:

   .. code-block:: cpp

      enum class SolverId {
          // ...
          MyPaper2026OneSided,
      };

3. **Add a case in the dispatch** in ``HomLib/robust/solver_dispatch.h``,
   inside ``dispatch_point_solver``:

   .. code-block:: cpp

      case SolverId::MyPaper2026OneSided:
          return func(MyPaper2026::SolverOneSided{});

4. **Add Python bindings** in ``python/src/main.cpp``:

   - Add a ``.value()`` to the ``SolverId`` enum binding.

   In ``python/src/homlib/__init__.py``:

   - Add a convenience class:

   .. code-block:: python

      class MyPaper2026OneSided(Solver):
          """One-sided solver by Author (Venue 2026), N-point."""
          def __init__(self):
              super().__init__(SolverId.MyPaper2026OneSided)

   - Add it to ``__all__``.

5. **Add tests** in ``python/tests/`` and/or ``tests/``.


Adding a non-minimal solver
----------------------------

Non-minimal solvers are solvers whose ``solve()`` natively handles any number
of points (not just the minimal sample). They are used for overdetermined
estimation inside RANSAC's local optimisation.

The steps are **identical** to adding a minimal solver. Once registered in the
dispatch, the solver is automatically available as a ``non_minimal_solver``
parameter in all backends:

.. code-block:: python

   # Use Wadenback's minimal solver with Nakano's overdetermined solver
   lomsac = homlib.LOMSAC(
       "Wadenback3DV2026OneSided",
       non_minimal_strategy="overdetermined",
       non_minimal_solver="NakanoICPR2025OneSided",
       inlier_threshold=0.07,
   )

The non-minimal solver's ``distortion_case`` must match the minimal solver's.
This is validated at construction time.

The type-erased function ``resolve_non_minimal_fn()`` in ``solver_dispatch.h``
wraps any registered solver's ``estimate()`` into a ``NonMinimalEstimateFn``,
so no additional wiring is needed.


One-sided-right solvers
-----------------------

If your solver is one-sided-left (distortion on the first image), HomLib can
automatically derive the right-sided variant using ``FlippedOneSidedSolver<T>``
in ``flipped_solver.h``. Add a second ``SolverId`` and dispatch case:

.. code-block:: cpp

   case SolverId::MyPaper2026OneSidedRight:
       return func(FlippedOneSidedSolver<MyPaper2026::SolverOneSided>{});

The flipped solver swaps the point arguments, inverts the homography, and
sets ``DISTORTION_CASE = ONE_SIDED_RIGHT``.
