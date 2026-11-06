"""
Synthetic benchmark of HomLib's point solvers, on problems from ``generate_problem_instance``.

For each solver and each distortion case it models (all solvers model NO_DISTORTION; the two-sided
solvers also model the one-sided cases), on minimal samples:
- without noise: the share of the problems where the solver misses the true solution (no solution
  within 1e-6, in homography error and |dk1| + |dk2|), and the median errors (numerical stability,
  as in Fig. 3 of Wadenbäck et al., 3DV 2026);
- with noise (``--noise``, in pixels for a focal length of 1000 px): the median error of the best
  solution (sensitivity to noise, as in their Fig. 4);
- the median time per call, from Python (it includes the overhead of the bindings).

``generate_problem_instance`` seeds itself with the time, so two runs don't solve the same problems:
``--compare`` reports statistical differences. ``--output`` and ``--compare`` otherwise work as in
hpatches_benchmark.py.

Usage
-----
    python synthetic_benchmark.py [--problems 2000] [--noise 0 0.5 1 2] [--solvers ...]
                                  [--output results.json]
    python synthetic_benchmark.py --compare before.json after.json
"""

import argparse
import json
import time

import numpy as np
from tqdm import tqdm

import homlib

FOCAL = 1000.0   # pixels: the noise is given in pixels for this focal length


def cases_of(name):
    """The distortion cases that a solver models"""
    C = homlib.DistortionCase
    if "TwoSidedEqual" in name:
        return [C.NO_DISTORTION, C.TWO_SIDED_EQUAL]
    if "TwoSided" in name:
        return [C.NO_DISTORTION, C.ONE_SIDED_LEFT, C.ONE_SIDED_RIGHT, C.TWO_SIDED_EQUAL, C.TWO_SIDED]
    if "OneSidedRight" in name:
        return [C.NO_DISTORTION, C.ONE_SIDED_RIGHT]
    if "OneSided" in name:
        return [C.NO_DISTORTION, C.ONE_SIDED_LEFT]
    return [C.NO_DISTORTION]


def sample_size(name):
    if "6Pt" in name:
        return 6
    if name in ("CaiPAMI2025", "DLT", "GuoPRL2019"):
        return 4
    return 5


def best_error(instance, solutions):
    """The error of the best solution: max(homography error, |dk1| + |dk2|), inf without one"""
    errors = [max(instance.hom_error(np.array(p.homography)),
                  2 * instance.dist_error(p.distortion_parameter, p.distortion_parameter2))
              for p in solutions]
    return min([e for e in errors if np.isfinite(e)], default=np.inf)


def point_solvers():
    return [s for s in homlib.SolverId.__members__.values()
            if homlib.Solver(s).feature_type == homlib.FeatureType.Points]


def main():
    parser = argparse.ArgumentParser(description="Synthetic benchmark of HomLib's point solvers")
    parser.add_argument("--problems", type=int, default=2000, help="Problems per solver, case and noise")
    parser.add_argument("--noise", type=float, nargs="*", default=[0.0, 0.5, 1.0, 2.0],
                        help="Noise levels in pixels (for a focal length of 1000 px)")
    parser.add_argument("--solvers", nargs="*", default=None, help="SolverId names (default: every point solver)")
    parser.add_argument("--output", type=str, default=None, help="Save the results to this JSON file")
    parser.add_argument("--compare", nargs=2, metavar=("A", "B"), default=None,
                        help="Print the differences B - A between two saved runs, and exit")
    args = parser.parse_args()
    if args.compare:
        compare(*args.compare)
        return

    solvers = ([homlib.SolverId.__members__[s] for s in args.solvers] if args.solvers
               else point_solvers())
    # results[solver][case][noise] = per problem: (best error, number of solutions, time)
    results = {}
    for s in tqdm(solvers, desc="Solvers"):
        solver = homlib.Solver(s)
        n = sample_size(s.name)
        results[s.name] = {}
        for case in cases_of(s.name):
            results[s.name][case.name] = {}
            for noise in args.noise:
                runs = []
                for _ in range(args.problems):
                    instance = homlib.generate_problem_instance(homlib.ProblemConfig(case, noise / FOCAL, n))
                    x, y = np.array(instance.x1).T, np.array(instance.x2).T
                    t0 = time.perf_counter()
                    solutions = solver.estimate(np.ascontiguousarray(x), np.ascontiguousarray(y))
                    elapsed = time.perf_counter() - t0
                    runs.append((best_error(instance, solutions), len(solutions), elapsed))
                results[s.name][case.name][str(noise)] = runs
    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f)
    print_results(results)


def summary(runs):
    err, n_sol, t = (np.array(v, dtype=float) for v in zip(*runs))
    return dict(missed=np.mean(~(err <= 1e-6)) * 100, median_error=np.median(err),
                solutions=np.mean(n_sol), time=np.median(t) * 1e6)


def print_results(results):
    print(f"\n{'Solver':<36s} {'Case':<16s} {'Noise':>5s} {'Missed':>8s} {'Median error':>13s} "
          f"{'Solutions':>9s} {'Time (us)':>9s}")
    for name, per_case in sorted(results.items()):
        for case, per_noise in per_case.items():
            for noise, runs in per_noise.items():
                r = summary(runs)
                missed = f"{r['missed']:7.2f}%" if float(noise) == 0 else ""
                print(f"{name:<36s} {case:<16s} {float(noise):5.1f} {missed:>8s} {r['median_error']:13.3g} "
                      f"{r['solutions']:9.2f} {r['time']:9.1f}")


def compare(path_a, path_b):
    """The differences between two runs, b - a: missed true solutions (noise-free) and the ratio of
    the median errors and times. The problems are the same in both runs only if the generator is
    seeded the same way, which it isn't (it seeds itself with the time): the differences are
    statistical, over --problems problems."""
    with open(path_a) as f:
        a = json.load(f)
    with open(path_b) as f:
        b = json.load(f)
    print(f"{'Solver':<36s} {'Case':<16s} {'Noise':>5s} {'Missed (pts)':>12s} {'Error ratio':>11s} {'Time':>6s}")
    for name in sorted(set(a) & set(b)):
        for case in a[name]:
            for noise in a[name].get(case, {}):
                if noise not in b[name].get(case, {}):
                    continue
                ra, rb = summary(a[name][case][noise]), summary(b[name][case][noise])
                missed = f"{rb['missed'] - ra['missed']:+11.2f}" if float(noise) == 0 else ""
                print(f"{name:<36s} {case:<16s} {float(noise):5.1f} {missed:>12s} "
                      f"{rb['median_error'] / ra['median_error']:11.2f} {rb['time'] / ra['time']:5.2f}x")


if __name__ == "__main__":
    main()
