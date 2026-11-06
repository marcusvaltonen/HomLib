"""
Synthetic benchmark of HomLib's point solvers, on problems from ``generate_problem_instance``.

Each solver solves minimal samples of the distortion cases it models (see ``SOLVERS``). For each
case and noise level, the script reports:
- how often the solver returns no solution at all;
- without noise, how often it misses the true solution: none of its solutions is within 1e-6 of
  it, in homography error and |dk1| + |dk2| (numerical stability, as in Fig. 3 of Wadenbäck et al.,
  3DV 2026);
- the median error of its best solution, over the problems with a solution (sensitivity to noise,
  as in their Fig. 4, with ``--noise`` in pixels for a focal length of 1000 px);
- the median time of a call, measured from Python, so it includes the overhead of the bindings.

``generate_problem_instance`` seeds itself with the time, so two runs don't solve the same
problems: ``--compare`` reports statistical differences, with their standard errors.

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

# For each point solver: its minimal sample size, and the distortion cases it models. In
# generate_problem_instance, ONE_SIDED_LEFT distorts the second image (distortion_parameter2), the
# side that the OneSided solvers estimate, and ONE_SIDED_RIGHT the first.
_C = homlib.DistortionCase
_ALL = [_C.NO_DISTORTION, _C.ONE_SIDED_LEFT, _C.ONE_SIDED_RIGHT, _C.TWO_SIDED_EQUAL, _C.TWO_SIDED]
SOLVERS = {
    "CaiPAMI2025": (4, [_C.NO_DISTORTION]),
    "DLT": (4, [_C.NO_DISTORTION]),
    "GuoPRL2019": (4, [_C.NO_DISTORTION]),
    "FitzgibbonCVPR2001OneSided": (5, [_C.NO_DISTORTION, _C.ONE_SIDED_LEFT]),
    "FitzgibbonCVPR2001TwoSidedEqual": (5, [_C.NO_DISTORTION, _C.TWO_SIDED_EQUAL]),
    "KukelovaCVPR2015TwoSided": (5, _ALL),
    "KukelovaCVPR2015TwoSided6Pt": (6, _ALL),
    "KukelovaCVPR2015TwoSidedEqual": (5, [_C.NO_DISTORTION, _C.TWO_SIDED_EQUAL]),
    "KukelovaCVPR2015TwoSidedEqual6Pt": (6, [_C.NO_DISTORTION, _C.TWO_SIDED_EQUAL]),
    "NakanoICPR2025OneSided": (5, [_C.NO_DISTORTION, _C.ONE_SIDED_LEFT]),
    "NakanoICPR2025OneSidedRight": (5, [_C.NO_DISTORTION, _C.ONE_SIDED_RIGHT]),
    "Wadenback3DV2026OneSided": (5, [_C.NO_DISTORTION, _C.ONE_SIDED_LEFT]),
    "Wadenback3DV2026TwoSided": (5, _ALL),
    "Wadenback3DV2026TwoSidedEqual": (5, [_C.NO_DISTORTION, _C.TWO_SIDED_EQUAL]),
}


def best_error(instance, solutions):
    """The error of the best solution, max(homography error, |dk1| + |dk2|); inf without one."""
    errors = [max(instance.hom_error(np.array(p.homography)),
                  2 * instance.dist_error(p.distortion_parameter, p.distortion_parameter2))
              for p in solutions]
    return min([e for e in errors if np.isfinite(e)], default=np.inf)


def main():
    parser = argparse.ArgumentParser(description="Synthetic benchmark of HomLib's point solvers")
    parser.add_argument("--problems", type=int, default=2000, help="Problems per solver, case and noise")
    parser.add_argument("--noise", type=float, nargs="*", default=[0.0, 0.5, 1.0, 2.0],
                        help="Noise levels in pixels, for a focal length of 1000 px")
    parser.add_argument("--solvers", nargs="*", choices=sorted(SOLVERS), default=sorted(SOLVERS))
    parser.add_argument("--output", type=str, default=None, help="Save the results to this JSON file")
    parser.add_argument("--compare", nargs=2, metavar=("A", "B"), default=None,
                        help="Print the differences B - A between two saved runs, and exit")
    args = parser.parse_args()
    if args.compare:
        compare(*args.compare)
        return

    # results[solver][case][noise] = per problem: (best error, number of solutions, time)
    results = {}
    for name in tqdm(args.solvers, desc="Solvers"):
        solver = homlib.Solver(homlib.SolverId.__members__[name])
        sample_size, cases = SOLVERS[name]
        results[name] = {}
        for case in cases:
            results[name][case.name] = {}
            for noise in args.noise:
                runs = []
                for _ in range(args.problems):
                    config = homlib.ProblemConfig(case, noise / FOCAL, sample_size)
                    instance = homlib.generate_problem_instance(config)
                    x, y = np.array(instance.x1).T, np.array(instance.x2).T
                    t0 = time.perf_counter()
                    solutions = solver.estimate(np.ascontiguousarray(x), np.ascontiguousarray(y))
                    elapsed = time.perf_counter() - t0
                    runs.append((best_error(instance, solutions), len(solutions), elapsed))
                results[name][case.name][str(noise)] = runs
    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f)
    print_results(results)


def summary(runs):
    err, n_sol, t = (np.array(v, dtype=float) for v in zip(*runs))
    solved = np.isfinite(err)
    return dict(n=len(err), none=np.mean(~solved), missed=np.mean(~(err <= 1e-6)),
                median_error=np.median(err[solved]) if solved.any() else np.inf,
                time=np.median(t) * 1e6)


def print_results(results):
    print(f"\n{'Solver':<34s} {'Case':<16s} {'Noise':>5s} {'No solution':>11s} {'Missed':>7s} "
          f"{'Median error':>12s} {'Time (us)':>9s}")
    for name, per_case in sorted(results.items()):
        for case, per_noise in per_case.items():
            for noise, runs in per_noise.items():
                r = summary(runs)
                missed = f"{r['missed'] * 100:6.2f}%" if float(noise) == 0 else ""
                print(f"{name:<34s} {case:<16s} {float(noise):5.1f} {r['none'] * 100:10.2f}% {missed:>7s} "
                      f"{r['median_error']:12.3g} {r['time']:9.1f}")


def difference(pa, pb, na, nb):
    """b - a for two rates, in percentage points, with its standard error."""
    se = np.sqrt(pa * (1 - pa) / na + pb * (1 - pb) / nb)
    return f"{(pb - pa) * 100:+6.2f} ±{se * 100:4.2f}"


def compare(path_a, path_b):
    """The differences between two runs, b - a: the rates of no solution and (without noise) of
    missed true solutions, in percentage points with their standard errors, and the ratios of the
    median errors and times."""
    with open(path_a) as f:
        a = json.load(f)
    with open(path_b) as f:
        b = json.load(f)
    print(f"{'Solver':<34s} {'Case':<16s} {'Noise':>5s} {'No solution':>13s} {'Missed':>13s} "
          f"{'Error':>6s} {'Time':>6s}")
    for name in sorted(set(a) & set(b)):
        for case, per_noise in a[name].items():
            for noise, runs in per_noise.items():
                if noise not in b[name].get(case, {}):
                    continue
                ra, rb = summary(runs), summary(b[name][case][noise])
                none = difference(ra["none"], rb["none"], ra["n"], rb["n"])
                missed = difference(ra["missed"], rb["missed"], ra["n"], rb["n"]) if float(noise) == 0 else ""
                finite = np.isfinite(ra["median_error"]) and np.isfinite(rb["median_error"])
                error = f"{rb['median_error'] / ra['median_error']:5.2f}x" if finite else "  n/a"
                print(f"{name:<34s} {case:<16s} {float(noise):5.1f} {none:>13s} {missed:>13s} {error:>6s} "
                      f"{rb['time'] / ra['time']:5.2f}x")


if __name__ == "__main__":
    main()
