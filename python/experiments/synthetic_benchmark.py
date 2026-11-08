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

``--outliers`` runs a second experiment, as in Fig. 5 of Wadenbäck et al.: LOMSAC on problems with
200 points, 0.5 px of noise and a share of outliers, with a 5 px inlier threshold. Each solver
solves problems of the distortion case it is made for (the last of its cases in ``SOLVERS``), so
the rows of different cases don't compare. For each solver and share, the script reports the share
of the true inliers that LOMSAC finds, and the median time of a run. Without ``--iterations``,
LOMSAC stops when it is confident or after 10000 iterations. ``--iterations N`` runs exactly N
iterations, so that every solver gets the same budget.

Every problem is seeded from ``--seed`` (default 0), the distortion case, the noise level (or the
share of outliers) and its index: the solvers with the same sample size solve the same problems,
and two runs with the same seed solve the same problems and give the same results on the same
platform (the C++ and NumPy random distributions differ between standard libraries), which helps
debugging. ``--compare`` reports the differences of the rates (no solution, missed, inliers found)
with their standard errors, and the ratios of the median errors and times, which have none; without
noise, the median errors are close to the machine precision, so their ratio is not given.
``--compare`` also warns about the settings that differ.

Usage
-----
    python synthetic_benchmark.py [--problems 2000] [--noise 0 0.5 1 2] [--solvers ...]
                                  [--output results.json]
    python synthetic_benchmark.py --outliers 0.2 0.4 0.6 0.8 [--iterations 100] [--problems 200]
                                  [--solvers ...] [--output outliers.json]
    python synthetic_benchmark.py --compare before.json after.json
"""

import argparse
import json
import time

import numpy as np
from tqdm import tqdm

import homlib
from hpatches_benchmark import warn_settings

FOCAL = 1000.0   # pixels: the noise is given in pixels for this focal length

# For each point solver: its minimal sample size, and the distortion cases it models, the last being
# the case it is made for (used by --outliers). In
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
_POINT_SOLVERS = {s.name for s in homlib.solver_ids(homlib.FeatureType.Points)}
assert set(SOLVERS) == _POINT_SOLVERS, ("SOLVERS must list the point solvers of homlib.SolverId: missing "
                                        f"{sorted(_POINT_SOLVERS - set(SOLVERS))}, "
                                        f"unknown {sorted(set(SOLVERS) - _POINT_SOLVERS)}")


def best_error(instance, solutions):
    """The error of the best solution, max(homography error, |dk1| + |dk2|); inf without one."""
    errors = [max(instance.hom_error(np.array(p.homography)),
                  2 * instance.dist_error(p.distortion_parameter, p.distortion_parameter2))
              for p in solutions]
    return min([e for e in errors if np.isfinite(e)], default=np.inf)


def main():
    parser = argparse.ArgumentParser(description="Synthetic benchmark of HomLib's point solvers")
    parser.add_argument("--problems", type=int, default=None,
                        help="Problems per solver, case and noise (default: 2000; with --outliers, 200)")
    parser.add_argument("--noise", type=float, nargs="+", default=None,
                        help="Noise levels in pixels, for a focal length of 1000 px (default: 0 0.5 1 2)")
    parser.add_argument("--solvers", nargs="+", choices=sorted(SOLVERS), default=sorted(SOLVERS))
    parser.add_argument("--outliers", type=float, nargs="+", default=None,
                        help="Shares of outliers: run LOMSAC on problems with outliers instead")
    parser.add_argument("--seed", type=int, default=0,
                        help="Seed of the problems (default: 0): same seed, same results on the same platform")
    parser.add_argument("--iterations", type=int, default=None,
                        help="With --outliers: run LOMSAC for exactly this many iterations")
    parser.add_argument("--output", type=str, default=None, help="Save the settings and the results to this JSON file")
    parser.add_argument("--compare", nargs=2, metavar=("A", "B"), default=None,
                        help="Print the differences B - A between two saved runs, and exit")
    args = parser.parse_args()
    if args.compare:
        compare(*args.compare)
        return

    if args.outliers is None and args.iterations is not None:
        parser.error("--iterations needs --outliers")
    if args.outliers is not None:
        if args.noise is not None:
            parser.error("--outliers runs with 0.5 px of noise: it doesn't take --noise")
        if not all(0 <= share < 1 for share in args.outliers):
            parser.error("the shares of outliers must be in [0, 1)")
        results = run_outliers(args.solvers, args.outliers, args.problems or 200, args.iterations, args.seed)
        if args.output:
            settings = dict(problems=args.problems or 200, noise=None, iterations=args.iterations, seed=args.seed,
                            outliers=args.outliers)
            with open(args.output, "w") as f:
                json.dump(dict(settings=settings, outliers=results), f)
        print_outliers(results)
        return

    # results[solver][case][noise] = per problem: (best error, number of solutions, time)
    results = {}
    problems = args.problems or 2000
    noises = args.noise or [0.0, 0.5, 1.0, 2.0]
    for name in tqdm(args.solvers, desc="Solvers"):
        solver = homlib.Solver(homlib.SolverId.__members__[name])
        sample_size, cases = SOLVERS[name]
        results[name] = {}
        for case in cases:
            results[name][case.name] = {}
            for noise in noises:
                runs = []
                for i in range(problems):
                    instance = make_problem(case, noise, sample_size,
                                            problem_seed(args.seed, int(case.value), round(noise * 1000), i))
                    x, y = np.array(instance.x1).T, np.array(instance.x2).T
                    t0 = time.perf_counter()
                    solutions = solver.estimate(np.ascontiguousarray(x), np.ascontiguousarray(y))
                    elapsed = time.perf_counter() - t0
                    runs.append((best_error(instance, solutions), len(solutions), elapsed))
                results[name][case.name][str(noise)] = runs
    if args.output:
        settings = dict(problems=problems, noise=noises, iterations=None, outliers=None, seed=args.seed)
        with open(args.output, "w") as f:
            json.dump(dict(settings=settings, results=results), f)
    print_results(results)


def load_results(path):
    """(settings, results) of a file saved with --output, results being {"outliers": results} with
    --outliers. The files saved before the settings have none: their settings are None."""
    with open(path) as f:
        data = json.load(f)
    settings = data.pop("settings", None)
    return settings, data.get("results", data)


def problem_seed(*key):
    """The seed of a problem, from (--seed, case, noise level or share of outliers, index): the value of
    the noise or share, not its position in the list, so that a level gets the same problems in every run"""
    return int(np.random.SeedSequence([int(k) for k in key]).generate_state(1)[0])


def make_problem(case, noise, n_points, seed):
    config = homlib.ProblemConfig(case, noise / FOCAL, n_points)
    config.seed = seed
    return homlib.generate_problem_instance(config)


def run_outliers(solvers, shares, problems, iterations, seed, n_points=200, noise=0.5, threshold=5.0):
    """results[solver][share] = per problem: (share of the true inliers found, time)"""
    results = {}
    options = {} if iterations is None else dict(min_num_iterations=iterations, max_num_iterations=iterations)
    for name in tqdm(solvers, desc="Solvers"):
        solver = homlib.Solver(homlib.SolverId.__members__[name])
        case = SOLVERS[name][1][-1]   # the case the solver is made for
        results[name] = {}
        for share in shares:
            runs = []
            for i in range(problems):
                ps = problem_seed(seed, int(case.value), round(share * 1000), i)
                instance = make_problem(case, noise, n_points, ps)
                rng = np.random.default_rng(ps)
                x, y = np.array(instance.x1).T, np.array(instance.x2).T
                n_out = int(round(share * n_points))
                # Outliers: random points in the bounding box of the second image's points
                y[:, :n_out] = rng.uniform(y.min(1, keepdims=True), y.max(1, keepdims=True), (2, n_out))
                lomsac = homlib.LOMSAC(solver, inlier_threshold=threshold / FOCAL, random_seed=ps, **options)
                t0 = time.perf_counter()
                _, stats = lomsac.estimate(np.ascontiguousarray(x), np.ascontiguousarray(y))
                elapsed = time.perf_counter() - t0
                found = np.sum(np.array(stats.inlier_indices, dtype=int) >= n_out)
                runs.append((found / (n_points - n_out), elapsed))
            results[name][str(share)] = runs
    return results


def print_outliers(results):
    print(f"\n{'Solver':<34s} {'Case':<16s} {'Outliers':>8s} {'Inliers found':>13s} {'Time (ms)':>9s}")
    for name, per_share in sorted(results.items()):
        case = SOLVERS[name][1][-1].name
        for share, runs in per_share.items():
            found, t = (np.array(v) for v in zip(*runs))
            print(f"{name:<34s} {case:<16s} {float(share) * 100:7.0f}% {np.mean(found) * 100:12.2f}% "
                  f"{np.median(t) * 1e3:9.2f}")


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
    median errors (with noise) and times, which have no standard error."""
    settings_a, a = load_results(path_a)
    settings_b, b = load_results(path_b)
    if ("outliers" in a) != ("outliers" in b):
        raise SystemExit("can't compare a run with --outliers and a run without it")
    warn_settings(settings_a, settings_b)
    if "outliers" in a:
        compare_outliers(a["outliers"], b["outliers"])
        return
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
                # Without noise, the median errors are close to the machine precision: no ratio
                finite = float(noise) > 0 and np.isfinite(ra["median_error"]) and np.isfinite(rb["median_error"])
                error = f"{rb['median_error'] / ra['median_error']:5.2f}x" if finite else "  n/a"
                print(f"{name:<34s} {case:<16s} {float(noise):5.1f} {none:>13s} {missed:>13s} {error:>6s} "
                      f"{rb['time'] / ra['time']:5.2f}x")


def compare_outliers(a, b):
    """b - a: the share of the true inliers found, in percentage points with its standard error, and
    the ratio of the median times"""
    print(f"{'Solver':<34s} {'Case':<16s} {'Outliers':>8s} {'Inliers found':>20s} {'Time':>7s}")
    for name in sorted(set(a) & set(b)):
        case = SOLVERS[name][1][-1].name
        for share in a[name]:
            if share not in b[name]:
                continue
            fa, ta = (np.array(v) for v in zip(*a[name][share]))
            fb, tb = (np.array(v) for v in zip(*b[name][share]))
            se = np.sqrt(fa.var(ddof=1) / len(fa) + fb.var(ddof=1) / len(fb)) * 100
            found = f"{(fb.mean() - fa.mean()) * 100:+7.2f} ±{se:5.2f} pts"
            print(f"{name:<34s} {case:<16s} {float(share) * 100:7.0f}% {found:>20s} "
                  f"{np.median(tb) / np.median(ta):6.3g}x")


if __name__ == "__main__":
    main()
