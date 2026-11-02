"""
Exhaustive SupeRANSAC configuration sweep on real image data.

Tests all combinations of sampler, scoring, local optimization, and final
optimization using multiprocessing with per-config timeouts.
"""

import random
import ctypes
import multiprocessing as mp
import time
import sys

import cv2
import numpy as np

import homlib

assert homlib._HAS_SUPERANSAC, "Requires SupeRANSAC build"

# ---------- Data setup (graf images) ----------

IMAGE1_PATH = 'img/grafA.png'
IMAGE2_PATH = 'img/grafB.png'
GT_HOMOGRAPHY_PATH = 'img/graf_model.txt'
DIST_COEFF_GT = -0.000003
SNN_THRESHOLD = 0.8

img1 = cv2.cvtColor(cv2.imread(IMAGE1_PATH), cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(cv2.imread(IMAGE2_PATH), cv2.COLOR_BGR2RGB)
K = np.eye(3); K[0, 2] = img2.shape[1] / 2.0; K[1, 2] = img2.shape[0] / 2.0
map_x, map_y = cv2.initUndistortRectifyMap(
    K, np.array([0, 0, 0, 0, 0, DIST_COEFF_GT, 0, 0]),
    np.eye(3), K, (img2.shape[1], img2.shape[0]), cv2.CV_32FC1)
img1 = cv2.remap(img1, map_x, map_y, cv2.INTER_LINEAR)
H_gt = np.linalg.inv(np.loadtxt(GT_HOMOGRAPHY_PATH))
H_gt_norm = H_gt / np.linalg.norm(H_gt)

w1, h1 = img1.shape[1], img1.shape[0]
w2, h2 = img2.shape[1], img2.shape[0]
T1 = np.eye(3); T1[0, 2] = -w1 / 2.0; T1[1, 2] = -h1 / 2.0
T2_inv = np.eye(3); T2_inv[0, 2] = w2 / 2.0; T2_inv[1, 2] = h2 / 2.0

# ---------- Feature extraction ----------

sift = cv2.SIFT_create(8000)
kps1, descs1 = sift.detectAndCompute(img1, None)
kps2, descs2 = sift.detectAndCompute(img2, None)
knn_matches = cv2.BFMatcher().knnMatch(descs1, descs2, k=2)

snn_ratios, tentatives = [], []
for m, n in knn_matches:
    if m.distance < SNN_THRESHOLD * n.distance:
        tentatives.append(m)
        snn_ratios.append(m.distance / n.distance)
sorted_indices = np.argsort(snn_ratios)
tentatives = list(np.array(tentatives)[sorted_indices])

src_pts = np.float32([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 2).T
dst_pts = np.float32([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 2).T
src_ori = np.float32([kps1[m.queryIdx].angle for m in tentatives])
dst_ori = np.float32([kps2[m.trainIdx].angle for m in tentatives])
ori = np.vstack((src_ori, dst_ori)) / 180.0 * np.pi

src_centered = src_pts - np.array([[w1 / 2.0], [h1 / 2.0]])
dst_centered = dst_pts - np.array([[w2 / 2.0], [h2 / 2.0]])

print(f"Matches: {len(tentatives)}")

# ---------- Helpers ----------

N_RUNS = 3
TIMEOUT_PER_CONFIG = 20  # seconds

# Enum maps
SAMP_MAP = {s: getattr(homlib.SuperansacSamplerType, s)
            for s in ["Uniform", "PROSAC", "NAPSAC", "ProgressiveNAPSAC",
                       "ImportanceSampler", "ARSampler"]}
SCORE_MAP = {s: getattr(homlib.SuperansacScoringType, s)
             for s in ["RANSAC", "MSAC", "MAGSAC", "MINPRAN", "ACRANSAC",
                        "GAU", "ML", "GRID"]}
LO_MAP = {s: getattr(homlib.SuperansacLocalOptType, s)
           for s in ["None_", "LSQ", "IRLS", "NestedRANSAC", "GCRANSAC",
                      "IteratedLMEDS", "CrossValidation"]}


def h_error(H_centered):
    H_pixel = T2_inv @ H_centered @ T1
    H_norm = H_pixel / np.linalg.norm(H_pixel)
    return np.linalg.norm(H_norm - H_gt_norm)


def run_one_config(args):
    """Run a single config in a worker process. Returns dict or error string."""
    solver_name, samp, score, lo, fo, thresh = args

    if solver_name == "Point":
        solver = homlib.Wadenback3DV2026OneSided()
        x, y, kw, invert = dst_centered, src_centered, {}, True
    else:
        solver = homlib.ValtonenOrnhagICPR2026OrientationSingleSided()
        x, y, kw, invert = src_centered, dst_centered, {"ori": ori}, False

    sa = homlib.SupeRANSAC(solver, inlier_threshold=thresh)
    sa.settings.sampler = SAMP_MAP[samp]
    sa.settings.scoring = SCORE_MAP[score]
    sa.settings.local_optimization = LO_MAP[lo]
    sa.settings.final_optimization = LO_MAP[fo]
    sa.settings.max_iterations = 3000
    sa.settings.min_iterations = 200

    errors, times = [], []
    for _ in range(N_RUNS):
        t0 = time.perf_counter()
        est, stats = sa.estimate(x, y, **kw)
        elapsed = time.perf_counter() - t0
        H = np.linalg.inv(est.homography) if invert else est.homography
        errors.append(h_error(H))
        times.append(elapsed * 1000)

    tag = f"{samp}|{score}|LO={lo}|F={fo}"
    return {
        "solver": solver_name,
        "tag": tag,
        "thresh": thresh,
        "H_err": float(np.median(errors)),
        "H_err_min": float(np.min(errors)),
        "ms": float(np.median(times)),
        "inliers": len(stats.inliers),
        "iters": int(stats.num_iterations),
    }


# ---------- LOMSAC baselines ----------

print("\n=== LOMSAC BASELINES ===")
for solver_name, solver, kw, x, y, invert in [
    ("Point", homlib.Wadenback3DV2026OneSided(), {}, dst_centered, src_centered, True),
    ("Orientation", homlib.ValtonenOrnhagICPR2026OrientationSingleSided(),
     {"ori": ori}, src_centered, dst_centered, False),
]:
    errors = []
    for _ in range(N_RUNS):
        opts = homlib.LORansacOptions()
        opts.squared_inlier_threshold = 9.0
        opts.final_least_squares = True
        opts.random_seed = random.randint(0, ctypes.c_uint32(-1).value)
        opts.lo_starting_iterations = 8
        opts.min_num_iterations = 70
        opts.max_num_iterations = 500
        lomsac = homlib.LOMSAC(solver, options=opts)
        est, stats = lomsac.estimate(x, y, **kw)
        H = np.linalg.inv(est.homography) if invert else est.homography
        errors.append(h_error(H))
    print(f"  {solver_name}: H_err median={np.median(errors):.4e}, "
          f"min={np.min(errors):.4e}, inliers={len(stats.inlier_indices)}")

# ---------- Build config grid ----------

SAMPLERS = list(SAMP_MAP.keys())
SCORINGS = list(SCORE_MAP.keys())
LOCAL_OPTS = list(LO_MAP.keys())
FINAL_OPTS = list(LO_MAP.keys())
THRESHOLDS = [10.0, 20.0]
SOLVER_NAMES = ["Point", "Orientation"]

configs = []
for solver_name in SOLVER_NAMES:
    for samp in SAMPLERS:
        for score in SCORINGS:
            for lo in LOCAL_OPTS:
                for fo in FINAL_OPTS:
                    for thresh in THRESHOLDS:
                        configs.append((solver_name, samp, score, lo, fo, thresh))

total = len(configs)
print(f"\n=== EXHAUSTIVE SupeRANSAC SWEEP ===")
print(f"Total configs: {total} "
      f"({len(SOLVER_NAMES)} solvers x {len(SAMPLERS)} samplers x "
      f"{len(SCORINGS)} scorings x {len(LOCAL_OPTS)} LO x "
      f"{len(FINAL_OPTS)} final x {len(THRESHOLDS)} thresholds)")
print(f"Workers: {mp.cpu_count()}, timeout: {TIMEOUT_PER_CONFIG}s per config")
print(flush=True)

# ---------- Run with process pool ----------

all_results = []
crashed = 0
timed_out = 0
crash_reasons = {}
done = 0

pool = mp.Pool(processes=min(mp.cpu_count(), 8))
async_results = []

for cfg in configs:
    ar = pool.apply_async(run_one_config, (cfg,))
    async_results.append((cfg, ar))

for cfg, ar in async_results:
    done += 1
    try:
        result = ar.get(timeout=TIMEOUT_PER_CONFIG)
        all_results.append(result)
    except mp.TimeoutError:
        timed_out += 1
        crash_reasons["TIMEOUT"] = crash_reasons.get("TIMEOUT", 0) + 1
    except Exception as e:
        crashed += 1
        reason = str(e).strip().split('\n')[-1][:80]
        crash_reasons[reason] = crash_reasons.get(reason, 0) + 1

    if done % 500 == 0:
        print(f"  [{done}/{total}] {len(all_results)} ok, "
              f"{crashed} crashed, {timed_out} timed out", flush=True)

pool.terminate()
pool.join()

print(f"\n  Done: {len(all_results)} ok, {crashed} crashed, {timed_out} timed out")

# ---------- Results per solver ----------

for solver_name in SOLVER_NAMES:
    results = sorted(
        [r for r in all_results if r["solver"] == solver_name],
        key=lambda r: r["H_err"])
    print(f"\n  TOP 20 for {solver_name}:")
    print(f"  {'Config':<55s} {'thr':>4s} {'H_err':>10s} {'min':>10s} "
          f"{'ms':>7s} {'inl':>5s} {'it':>5s}")
    print(f"  {'-' * 100}")
    for r in results[:20]:
        print(f"  {r['tag']:<55s} {r['thresh']:>4.0f} {r['H_err']:>10.2e} "
              f"{r['H_err_min']:>10.2e} {r['ms']:>7.1f} {r['inliers']:>5d} "
              f"{r['iters']:>5d}")

# ---------- Overall summary ----------

print(f"\n\n{'=' * 115}")
print(f"OVERALL TOP 30 (all solvers)")
print(f"{'=' * 115}")
all_results.sort(key=lambda r: r["H_err"])
print(f"{'Solver':<14s} {'Config':<55s} {'thr':>4s} {'H_err':>10s} "
      f"{'min':>10s} {'ms':>7s} {'inl':>5s}")
print("-" * 115)
for r in all_results[:30]:
    print(f"{r['solver']:<14s} {r['tag']:<55s} {r['thresh']:>4.0f} "
          f"{r['H_err']:>10.2e} {r['H_err_min']:>10.2e} {r['ms']:>7.1f} "
          f"{r['inliers']:>5d}")

# ---------- Best per axis ----------

print(f"\n\n=== BEST PER AXIS (threshold=10) ===")
t10 = [r for r in all_results if r["thresh"] == 10.0]
if t10:
    for axis_name, axis_key in [("sampler", 0), ("scoring", 1), ("LO", 2), ("final", 3)]:
        print(f"\n  Best {axis_name}:")
        by_val = {}
        for r in t10:
            parts = r["tag"].split("|")
            val = parts[axis_key].replace("LO=", "").replace("F=", "")
            by_val.setdefault(val, []).append(r["H_err"])
        for val, errs in sorted(by_val.items(), key=lambda x: np.median(x[1])):
            print(f"    {val:<20s} median={np.median(errs):.4e} "
                  f"best={np.min(errs):.4e} (n={len(errs)})")

print(f"\n\nTotal: {total} configs, "
      f"{len(all_results)} succeeded, {crashed} crashed, {timed_out} timed out")
if crash_reasons:
    print("\nFailure reasons:")
    for reason, count in sorted(crash_reasons.items(), key=lambda x: -x[1]):
        print(f"  [{count:>4d}x] {reason}")
