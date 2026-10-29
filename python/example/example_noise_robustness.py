r"""
Noise robustness
================

Evaluate how solver accuracy degrades with increasing point noise for
the two-sided (independent) distortion model. We simulate a camera with
focal length *f* = 1000 px and express noise in pixel standard deviation,
normalizing the points before passing them to the solvers.

The solvers compared are by Wadenbäck et al. [1]_ and Kukelova et al. [2]_.

References
----------
.. [1] M. Wadenbäck, M. Valtonen Örnhag, and J. Edstedt, "Radially Distorted
   Homographies, Revisited," in Proc. 3DV, 2026.
.. [2] Z. Kukelova, J. Heller, M. Bujnak, and T. Pajdla, "Radial Distortion
   Homography," in Proc. CVPR, 2015.
"""

import numpy as np
import matplotlib.pyplot as plt

import homlib

###############################################################################
# Setup
# -----
# Define focal length, noise sweep, and the two-sided equal solvers.

focal_length = 1000.0
noise_pixels = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
n_trials = 100

solver_specs = {
    "Wadenbäck (3DV 2026)": lambda x, y: homlib.estimate_wadenback_3dv_2026_two_sided(x, y, True),
    "Kukelova (CVPR 2015)": lambda x, y: homlib.estimate_kukelova_cvpr_2015_two_sided(x, y, True),
}
colors = ["#4c72b0", "#dd8452"]

###############################################################################
# Sweep noise levels
# ------------------
# For each noise level (in pixels) we convert to normalised image
# coordinates by dividing by the focal length.

all_errors = {name: [] for name in solver_specs}

for noise_px in noise_pixels:
    noise_norm = noise_px / focal_length

    per_solver = {name: [] for name in solver_specs}
    for _ in range(n_trials):
        config = homlib.ProblemConfig(
            homlib.DistortionCase.TWO_SIDED, noise_norm, 5
        )
        instance = homlib.generate_problem_instance(config)

        x = np.array(instance.x1).T
        y = np.array(instance.x2).T

        for name, fn in solver_specs.items():
            poses = fn(x, y)
            best_err = float("nan")
            for p in poses:
                err = instance.hom_error(p.homography)
                if not np.isfinite(best_err) or err < best_err:
                    best_err = err
            per_solver[name].append(best_err)

    for name in solver_specs:
        all_errors[name].append(np.nanmedian(per_solver[name]))

###############################################################################
# Print results

print(f"{'Noise (px)':>12s}", end="")
for name in solver_specs:
    print(f"   {name:>25s}", end="")
print()
print("-" * (12 + 28 * len(solver_specs)))
for i, noise_px in enumerate(noise_pixels):
    print(f"{noise_px:12.1f}", end="")
    for name in solver_specs:
        print(f"   {all_errors[name][i]:25.4e}", end="")
    print()

###############################################################################
# Homography error vs. noise
# ---------------------------
# Median homography error as a function of pixel noise.

plt.figure(figsize=(7, 4))
for (name, errs), c in zip(all_errors.items(), colors):
    plt.semilogy(noise_pixels, errs, "o-", label=name, color=c)
plt.xlabel("Point noise std (pixels, f = 1000)")
plt.ylabel("Median homography error")
plt.title("Noise robustness — two-sided distortion")
plt.legend(fontsize=9)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
