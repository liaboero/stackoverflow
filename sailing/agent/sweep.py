"""
Compare variants of the planner (parameter overrides on top of trained_params.json) on a development set of 48
episodes (seeds 90000+ and 91000+, disjoint from training and from the test set of evaluate.py). The comparison is
paired, episode by episode, against the first variant.

    python agent/sweep.py '{}' '{"gust_weight": 1}' '{"gust_weight": 1, "gust_decay": 20}'
"""
import json, os, sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import run_episode
from forecast_agent import ForecastPlanner

BASE = {k: v for k, v in json.load(open(os.path.join(HERE, "trained_params.json"))).items() if k != "train_score"}
DEV = [("train", 90000 + i) for i in range(24)] + [("random", 91000 + i) for i in range(24)]


def job(args):
    k, over, regime, seed = args
    return k, regime, seed, run_episode(ForecastPlanner({**BASE, **over}), seed, regime=regime)["score"]


if __name__ == "__main__":
    variants = [json.loads(a) for a in sys.argv[1:]] or [{}]
    with Pool(os.cpu_count()) as p:
        res = p.map(job, [(k, v, r, s) for k, v in enumerate(variants) for r, s in DEV], chunksize=1)
    sc = {(k, r, s): x for k, r, s, x in res}
    ref = np.array([sc[(0, r, s)] for r, s in DEV])
    for k, v in enumerate(variants):
        x = np.array([sc[(k, r, s)] for r, s in DEV])
        d = x - ref
        print(f"{json.dumps(v):60s} {x.mean():7.2f}   {d.mean():+.2f} ± {d.std(ddof=1) / np.sqrt(len(d)):.2f}", flush=True)
