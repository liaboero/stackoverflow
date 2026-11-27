"""
Large validation set for decisions about unknown maps: 400 episodes never used anywhere else, 300 on random regimes
(other wind settings and island layouts) and 100 on the training regime. Results are cached per parameter file, and
every comparison is paired (same episodes) against the first file given.

    python agent/validate.py agent/best_params.json agent/best_params_retrained.json
"""
import hashlib, json, os, sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import run_episode
from forecast_agent import ForecastPlanner

VAL = [("random", 120000 + i) for i in range(300)] + [("train", 110000 + i) for i in range(100)]
CACHE = os.path.join(HERE, "validation_cache")


def load(path):
    p = json.load(open(path))
    p.pop("train_score", None)
    return p


def job(args):
    params, regime, seed = args
    return regime, seed, run_episode(ForecastPlanner(params), seed, regime=regime)["score"]


def scores(params, pool):
    key = hashlib.sha1(json.dumps(params, sort_keys=True).encode()).hexdigest()[:16]
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, key + ".json")
    if os.path.exists(path):
        return {tuple(k.split(":", 1)[0:1]) + (int(k.split(":")[1]),): v for k, v in json.load(open(path)).items()}
    res = pool.map(job, [(params, r, s) for r, s in VAL], chunksize=1)
    out = {(r, s): x for r, s, x in res}
    json.dump({f"{r}:{s}": x for (r, s), x in out.items()}, open(path, "w"))
    return out


def summary(sc, ref=None):
    x = np.array([sc[e] for e in VAL])
    rnd, trn = x[:300].mean(), x[300:].mean()
    line = f"random {rnd:6.2f}  train {trn:6.2f}  all {x.mean():6.2f} ± {x.std(ddof=1) / np.sqrt(len(x)):.2f}"
    if ref is not None:
        d = x - np.array([ref[e] for e in VAL])
        line += f"   paired diff {d.mean():+.2f} ± {d.std(ddof=1) / np.sqrt(len(d)):.2f}"
    return line


if __name__ == "__main__":
    with Pool(os.cpu_count()) as pool:
        ref = None
        for f in sys.argv[1:]:
            sc = scores(load(f), pool)
            print(f"{os.path.basename(f):32s} {summary(sc, ref)}", flush=True)
            ref = ref or sc
