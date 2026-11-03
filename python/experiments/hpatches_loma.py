"""
HPatches benchmark with LoMa feature matching and HomLib robust estimation.

Compares LOMSAC and SupeRANSAC backends on HPatches image pairs using
LoMa (Nordström et al., ECCV 2026) for local feature matching.

Requirements
------------
pip install loma opencv-python matplotlib tqdm

HPatches dataset is downloaded automatically on first run.

Usage
-----
    python hpatches_loma.py [--data_dir /path/to/hpatches] [--num_keypoints 2048]
"""

import argparse
import time
from collections import defaultdict
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import homlib

HPATCHES_URL = "http://icvl.ee.ic.ac.uk/vgg/hpatches/sequences/hpatches-sequences-release.tar.gz"
HPATCHES_DIRNAME = "hpatches-sequences-release"


def download_hpatches(data_dir: Path) -> Path:
    """Download and extract HPatches if not already present."""
    seq_dir = data_dir / HPATCHES_DIRNAME
    if seq_dir.exists():
        return seq_dir

    import tarfile
    import urllib.request

    archive = data_dir / "hpatches-sequences-release.tar.gz"
    if not archive.exists():
        print(f"Downloading HPatches to {archive} ...")
        data_dir.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(HPATCHES_URL, archive)

    print("Extracting ...")
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(data_dir)

    return seq_dir


def load_hpatches_pairs(seq_dir: Path, max_sequences: int | None = None):
    """Yield (img_ref_path, img_i_path, H_gt) tuples from HPatches."""
    sequences = sorted(p for p in seq_dir.iterdir() if p.is_dir())
    if max_sequences is not None:
        sequences = sequences[:max_sequences]

    for seq in sequences:
        ref_path = seq / "1.ppm"
        if not ref_path.exists():
            continue
        for i in range(2, 7):
            img_path = seq / f"{i}.ppm"
            h_path = seq / f"H_1_{i}"
            if not img_path.exists() or not h_path.exists():
                continue
            H_gt = np.loadtxt(h_path).reshape(3, 3)
            yield seq.name, i, str(ref_path), str(img_path), H_gt


def center_keypoints(kpts, width, height):
    """Shift keypoints so (0,0) is at the image centre."""
    shift = np.array([[width / 2.0], [height / 2.0]])
    return kpts - shift


def corner_error(H_est, H_gt, width, height):
    """Mean corner reprojection error in pixels."""
    corners = np.array([[0, 0, 1], [width, 0, 1],
                        [width, height, 1], [0, height, 1]], dtype=np.float64).T
    proj_est = H_est @ corners
    proj_est /= proj_est[2:3]
    proj_gt = H_gt @ corners
    proj_gt /= proj_gt[2:3]
    return np.mean(np.linalg.norm(proj_est[:2] - proj_gt[:2], axis=0))


def run_lomsac(solver, x, y, threshold):
    lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=threshold ** 2)
    t0 = time.perf_counter()
    estimate, stats = lomsac.estimate(x, y)
    elapsed = time.perf_counter() - t0
    return estimate, stats, elapsed


def run_superansac(solver, x, y, threshold):
    superansac = homlib.SupeRANSAC(solver, inlier_threshold=threshold)
    t0 = time.perf_counter()
    estimate, stats = superansac.estimate(x, y)
    elapsed = time.perf_counter() - t0
    return estimate, stats, elapsed


def main():
    parser = argparse.ArgumentParser(description="HPatches + LoMa + HomLib benchmark")
    parser.add_argument("--data_dir", type=str, default="data",
                        help="Directory for HPatches data")
    parser.add_argument("--num_keypoints", type=int, default=2048,
                        help="Number of LoMa keypoints")
    parser.add_argument("--threshold", type=float, default=3.0,
                        help="Inlier threshold in pixels")
    parser.add_argument("--max_sequences", type=int, default=None,
                        help="Limit number of sequences (for quick testing)")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    seq_dir = download_hpatches(data_dir)

    # Load LoMa
    from loma import LoMa, LoMaB
    model = LoMa(LoMaB(num_keypoints=args.num_keypoints))

    has_superansac = homlib._HAS_SUPERANSAC

    solvers = {
        "Wadenback (one-sided)": homlib.Wadenback3DV2026OneSided(),
        "Fitzgibbon (one-sided)": homlib.FitzgibbonCVPR2001OneSided(),
        "Nakano (one-sided)": homlib.NakanoICPR2025OneSided(),
    }

    backends = ["LOMSAC"]
    if has_superansac:
        backends.append("SupeRANSAC")
    else:
        print("Note: SupeRANSAC not available, benchmarking LOMSAC only.")

    results = defaultdict(lambda: defaultdict(list))

    pairs = list(load_hpatches_pairs(seq_dir, args.max_sequences))
    print(f"Running on {len(pairs)} image pairs, {len(solvers)} solvers, "
          f"{len(backends)} backend(s)\n")

    for seq_name, pair_idx, ref_path, img_path, H_gt in tqdm(pairs, desc="Pairs"):
        # Get matches from LoMa
        try:
            kptsA, kptsB = model.match(ref_path, img_path)
        except Exception as e:
            print(f"  LoMa failed on {seq_name}/{pair_idx}: {e}")
            continue

        if len(kptsA) < 10:
            continue

        # Read image sizes for centering
        img_ref = cv2.imread(ref_path)
        img_tgt = cv2.imread(img_path)
        h1, w1 = img_ref.shape[:2]
        h2, w2 = img_tgt.shape[:2]

        # Center keypoints (HomLib expects origin at image centre)
        x = center_keypoints(kptsA.T, w1, h1)
        y = center_keypoints(kptsB.T, w2, h2)

        # Centering transforms for denormalizing the homography
        T1 = np.eye(3)
        T1[0, 2], T1[1, 2] = -w1 / 2.0, -h1 / 2.0
        T2_inv = np.eye(3)
        T2_inv[0, 2], T2_inv[1, 2] = w2 / 2.0, h2 / 2.0

        for solver_name, solver in solvers.items():
            # LOMSAC
            try:
                est_l, stats_l, t_l = run_lomsac(solver, x, y, args.threshold)
                H_l = T2_inv @ est_l.homography @ T1
                err_l = corner_error(H_l, H_gt, w1, h1)
            except Exception:
                err_l, t_l = float("nan"), float("nan")

            key_l = f"{solver_name}/LOMSAC"
            results[key_l]["errors"].append(err_l)
            results[key_l]["times"].append(t_l)

            # SupeRANSAC
            if has_superansac:
                try:
                    est_s, stats_s, t_s = run_superansac(
                        solver, x, y, args.threshold)
                    H_s = T2_inv @ est_s.homography @ T1
                    err_s = corner_error(H_s, H_gt, w1, h1)
                except Exception:
                    err_s, t_s = float("nan"), float("nan")

                key_s = f"{solver_name}/SupeRANSAC"
                results[key_s]["errors"].append(err_s)
                results[key_s]["times"].append(t_s)

    # -------------------------------------------------------------------------
    # Print summary
    # -------------------------------------------------------------------------
    auc_thresholds = [1, 3, 5, 10]

    print(f"\n{'Method':<45s}", end="")
    for th in auc_thresholds:
        print(f"  AUC@{th}px", end="")
    print(f"  {'Med. err (px)':>13s}  {'Med. time (ms)':>14s}")
    print("-" * (45 + 10 * len(auc_thresholds) + 30))

    for key in sorted(results.keys()):
        errors = np.array(results[key]["errors"])
        times = np.array(results[key]["times"])
        valid = np.isfinite(errors)
        errors_valid = errors[valid]

        print(f"{key:<45s}", end="")
        for th in auc_thresholds:
            auc = np.mean(np.minimum(errors_valid, th) / th) * 100 if len(errors_valid) else 0
            print(f"  {auc:7.2f}  ", end="")
        med_err = np.nanmedian(errors) if len(errors) else float("nan")
        med_time = np.nanmedian(times) * 1000 if len(times) else float("nan")
        print(f"  {med_err:13.2f}  {med_time:14.3f}")

    # -------------------------------------------------------------------------
    # Cumulative error plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(9, 5))
    colors = plt.cm.tab10.colors
    for i, key in enumerate(sorted(results.keys())):
        errors = np.sort(np.array(results[key]["errors"]))
        errors = errors[np.isfinite(errors)]
        if len(errors) == 0:
            continue
        cdf = np.arange(1, len(errors) + 1) / len(errors)
        plt.plot(errors, cdf, label=key, color=colors[i % len(colors)])

    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel("Corner error (px)")
    plt.ylabel("Fraction of pairs")
    plt.title("HPatches — LoMa matches, HomLib robust estimation")
    plt.legend(fontsize=7, loc="lower right")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("hpatches_loma_results.png", dpi=150)
    print("\nPlot saved to hpatches_loma_results.png")

    # -------------------------------------------------------------------------
    # Timing comparison bar chart
    # -------------------------------------------------------------------------
    if has_superansac:
        solver_names = list(solvers.keys())
        lomsac_times = [np.nanmedian(results[f"{s}/LOMSAC"]["times"]) * 1000
                        for s in solver_names]
        superansac_times = [np.nanmedian(results[f"{s}/SupeRANSAC"]["times"]) * 1000
                           for s in solver_names]

        x_pos = np.arange(len(solver_names))
        width = 0.35

        fig, ax = plt.subplots(figsize=(9, 4))
        ax.bar(x_pos - width / 2, lomsac_times, width, label="LOMSAC")
        ax.bar(x_pos + width / 2, superansac_times, width, label="SupeRANSAC")
        ax.set_ylabel("Median time (ms)")
        ax.set_title("LOMSAC vs SupeRANSAC — per-pair estimation time")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(solver_names, fontsize=8)
        ax.legend()
        ax.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig("hpatches_loma_timing.png", dpi=150)
        print("Timing plot saved to hpatches_loma_timing.png")

    plt.show()


if __name__ == "__main__":
    main()
