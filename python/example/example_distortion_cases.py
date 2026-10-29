r"""
Distortion models
=================

HomLib supports several radial distortion models. This example generates
a synthetic problem for each distortion case and runs the appropriate
solver, demonstrating the different models available.

Distortion conventions:

- **One-sided**: distortion on one image only
- **Two-sided equal**: same distortion on both images
- **Two-sided**: independent distortion on each image

All solvers used here are from Wadenbäck et al. [1]_.

References
----------
.. [1] M. Wadenbäck, M. Valtonen Örnhag, and J. Edstedt, "Radially Distorted
   Homographies, Revisited," in Proc. 3DV, 2026.
"""

import numpy as np

import homlib

###############################################################################
# Define a solver for each distortion case
#
# "Left" and "right" refer to the side of the homography equation
# ``undistort(y) = H * x``. ``ONE_SIDED_LEFT`` puts distortion on x2 (y),
# matching the one-sided solver which undistorts its second argument.

cases = [
    (
        "One-sided",
        homlib.DistortionCase.ONE_SIDED_LEFT,
        lambda x, y: homlib.estimate_wadenback_3dv_2026_one_sided(x, y, True),
    ),
    (
        "Two-sided equal",
        homlib.DistortionCase.TWO_SIDED_EQUAL,
        lambda x, y: homlib.estimate_wadenback_3dv_2026_two_sided_equal(x, y, True),
    ),
    (
        "Two-sided",
        homlib.DistortionCase.TWO_SIDED,
        lambda x, y: homlib.estimate_wadenback_3dv_2026_two_sided(x, y, True),
    ),
]

###############################################################################
# Run each case

print(f"{'Case':<22s} {'H error':>10s} {'k error':>10s} {'k1 est':>10s} {'k2 est':>10s} {'k1 gt':>10s} {'k2 gt':>10s}")
print("-" * 88)

for name, dist_case, solver_fn in cases:
    config = homlib.ProblemConfig(dist_case, 0.0, 5)
    instance = homlib.generate_problem_instance(config)

    x = np.array(instance.x1).T
    y = np.array(instance.x2).T

    poses = solver_fn(x, y)

    gt = instance.posedata

    best_err = float("inf")
    best_pose = None
    for p in poses:
        err = instance.hom_error(p.homography)
        if err < best_err:
            best_err = err
            best_pose = p

    dist_err = instance.dist_error(best_pose.distortion_parameter,
                                   best_pose.distortion_parameter2)
    print(
        f"{name:<22s} {best_err:10.2e} {dist_err:10.2e}"
        f" {best_pose.distortion_parameter:10.4f} {best_pose.distortion_parameter2:10.4f}"
        f" {gt.distortion_parameter:10.4f} {gt.distortion_parameter2:10.4f}"
    )
