r"""
Range check on the distortion
=============================

Minimal solvers that find the distortion as the roots of a polynomial return many
solutions, most with implausible distortion parameters. A range check on the
distortion keeps only the solutions with both parameters in ``[k_min, k_max]``. The solvers
apply it to each root, before computing its homography: RANSAC has fewer hypotheses to score,
and the solvers skip the homographies of the rejected roots.

The range is in HomLib's normalized coordinates: the robust estimators (LOMSAC and the other
backends) divide both images' points of each sample by
:math:`f_0 = \mathrm{mean}(|x|) / \sqrt{2}`, computed from the first image's points (measured
from the distortion center), whatever the units of the input. They
apply ``[-1, 0.1]`` by default; :meth:`~homlib.Solver.set_k_range` changes it, and
``set_k_range()`` without bounds disables it. The direct :meth:`~homlib.Solver.estimate` calls
the solver on the points as given, without a range unless one is set.

The default is safe for every lens that the division model can represent: in these coordinates,
for matches spread over the image, :math:`k` is about the relative distortion at the image corners
divided by 7, so ``[-1, 0.1]`` only rejects more than about 70% of pincushion distortion there.
Smartphones and action cameras, whose distortion López-Antequera et al. [1]_ measured on more than
300 cameras, and long lenses, with a few percent of pincushion distortion, are far inside (see the
solver guide, "Range check on the distortion").

The check only selects the hypotheses: the final refinement isn't constrained.

References
----------
.. [1] M. López-Antequera, R. Marí, P. Gargallo, Y. Kuang, J. Gonzalez-Jimenez, and
   G. Haro, "Deep Single Image Camera Calibration with Radial Distortion," in Proc. IEEE/CVF
   Conference on Computer Vision and Pattern Recognition (CVPR), 2019, pp. 11809-11817.
   https://research.mapillary.com/publication/cvpr19d
.. [2] Z. Kukelova, J. Heller, M. Bujnak, and T. Pajdla, "Radial Distortion Homography," in
   Proc. IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2015.
"""

import time

import numpy as np

import homlib

rng = np.random.default_rng()   # random problems: the numbers change from run to run

###############################################################################
# Fewer solutions per minimal sample
# ----------------------------------
#
# Kukelova et al.'s 5-point solver [2]_, called directly on noise-free problems with both
# images distorted, without and with the range ``[-1, 0.1]``.

solver = homlib.KukelovaCVPR2015TwoSided()
n_all = n_kept = n_true_kept = 0
n_problems = 200
for _ in range(n_problems):
    instance = homlib.generate_problem_instance(homlib.ProblemConfig(homlib.DistortionCase.TWO_SIDED, 0.0, 5))
    x = np.array(instance.x1).T
    y = np.array(instance.x2).T
    k1, k2 = instance.posedata.distortion_parameter, instance.posedata.distortion_parameter2

    solver.reset_k_range()   # estimate(): no range by default
    n_all += len(solver.estimate(x, y))
    solver.set_k_range(-1.0, 0.1)
    kept = solver.estimate(x, y)
    n_kept += len(kept)
    n_true_kept += any(abs(p.distortion_parameter - k1) < 1e-6 and abs(p.distortion_parameter2 - k2) < 1e-6
                       for p in kept)

print(f"Solutions per sample: {n_all / n_problems:.2f} without the range, {n_kept / n_problems:.2f} with it")
print(f"True solution kept: {n_true_kept} of {n_problems} problems")

###############################################################################
# Faster RANSAC
# -------------
#
# LOMSAC on 500 noisy matches, 40% of them outliers, with its default range and without a
# range (times: median of 5 runs). The points are in pixels, for a focal length of 1000 px:
# the range doesn't depend on their units.

instance = homlib.generate_problem_instance(homlib.ProblemConfig(homlib.DistortionCase.TWO_SIDED, 1e-3, 500))
f = 1000.0
x = f * np.array(instance.x1).T
y = f * np.array(instance.x2).T
n_outliers = 200
x[:, :n_outliers] = rng.uniform(x.min(), x.max(), (2, n_outliers))
y[:, :n_outliers] = rng.uniform(y.min(), y.max(), (2, n_outliers))

for label, use_range in [("with the default range", True), ("without a range", False)]:
    solver = homlib.KukelovaCVPR2015TwoSided()
    if not use_range:
        solver.set_k_range()   # no bounds: no check
    lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=5.0 ** 2)   # 5 px, 5 noise sigmas
    times = []
    for _ in range(5):
        start = time.perf_counter()
        estimate, stats = lomsac.estimate(x, y)
        times.append(time.perf_counter() - start)
    elapsed = np.median(times)
    correct = np.sum(np.array(stats.inlier_indices) >= n_outliers)
    print(f"LOMSAC {label}: {elapsed * 1e3:.1f} ms, {len(stats.inlier_indices)} inliers "
          f"({correct} of the {x.shape[1] - n_outliers} true ones), "
          f"k1, k2 (in the generator's units, i.e. pixels / 1000) = {estimate.distortion_parameter * f**2:.4f}, "
          f"{estimate.distortion_parameter2 * f**2:.4f} (true {instance.posedata.distortion_parameter:.4f}, "
          f"{instance.posedata.distortion_parameter2:.4f})")
