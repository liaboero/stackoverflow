"""Targeted test of one parameter, front_weight, on episodes WITH a turning front and episodes without one."""
import json, os, sys
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import SailingEnv, run_episode
from forecast_agent import ForecastPlanner

BASE = {k: v for k, v in json.load(open(os.path.join(HERE, "trained_params.json"))).items() if k != "train_score"}


def kind(regime, seed):
    env = SailingEnv(regime); env.reset(seed=seed)
    f = env.wind.front
    return "turning" if f and f["kind"] == "turning" else ("still" if f else "none")


def job(args):
    fw, regime, seed = args
    return fw, regime, seed, run_episode(ForecastPlanner({**BASE, "front_weight": fw}), seed, regime=regime)["score"]


if __name__ == "__main__":
    # seeds disjoint from training (1000-1131, 10000-11015) and from the test set (50000+, 60000+, scenarios)
    pool_eps = [("train", 70000 + i) for i in range(80)] + [("random", 80000 + i) for i in range(80)]
    kinds = {e: kind(*e) for e in pool_eps}
    turning = [e for e in pool_eps if kinds[e] == "turning"][:28]
    other = [e for e in pool_eps if kinds[e] != "turning"][:28]
    weights = [1.0, 1.5, 2.0, 2.5, 3.0]
    with Pool(os.cpu_count()) as p:
        res = p.map(job, [(w, r, s) for w in weights for r, s in turning + other], chunksize=1)
    sc = {(w, r, s): x for w, r, s, x in res}
    print(f"{'front_weight':>12s} {'turning front':>14s} {'others':>8s} {'all':>8s}   (n = {len(turning)} + {len(other)})")
    for w in weights:
        t = np.mean([sc[(w, *e)] for e in turning]); o = np.mean([sc[(w, *e)] for e in other])
        print(f"{w:12.1f} {t:14.2f} {o:8.2f} {(t + o) / 2:8.2f}")
    ref = np.array([sc[(1.5, *e)] for e in turning])
    for w in weights:
        d = np.array([sc[(w, *e)] for e in turning]) - ref
        print(f"turning, {w} vs 1.5: {d.mean():+.2f} ± {d.std(ddof=1) / np.sqrt(len(d)):.2f}")
