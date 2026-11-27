"""
Training of the ForecastPlanner parameters by random search with a shrinking step, on a fixed set of training
episodes (common random numbers: every candidate sails the same episodes). The test episodes (evaluate.py) never
appear here.

    python agent/train.py --gens 12 --pop 8 --out agent/params.json
"""
import argparse, json, os, sys, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from sailing_challenge import run_episode
from forecast_agent import ForecastPlanner, DEFAULT_PARAMS

# training episodes: half on the training regime, half on random regimes (seeds disjoint from evaluate.py)
TRAIN_EPISODES = [("train", 1000 + i) for i in range(16)] + [("random", 10000 + i) for i in range(16)]
if os.environ.get("SAILING_TRAIN_EPISODES") == "64":   # longer training: twice as many episodes
    TRAIN_EPISODES += [("train", 1100 + i) for i in range(16)] + [("random", 11000 + i) for i in range(16)]

SPACE = {  # name: (low, high, integer)
    "osc_weight": (0.0, 1.5, False),
    "osc_min_steps": (4, 30, True),
    "trend_weight": (0.0, 1.5, False),
    "front_weight": (0.0, 4.0, False),   # was (0, 1.5): the first two trainings ended on the bound
    "front_threshold": (5.0, 60.0, False),
    "front_prior_T": (8.0, 80.0, False),
    "lag_frac": (0.0, 1.0, False),
}
FIXED = {"horizon": 5, "beam": 16}

# --space v2: the parameters of the agent with the gust forecast and the Bayesian oscillation forecast
SPACE_V2 = {
    "osc_weight": (0.0, 1.5, False),
    "osc_sigma": (1.0, 12.0, False),
    "front_weight": (0.0, 2.0, False),
    "front_threshold": (5.0, 60.0, False),
    "front_prior_T": (8.0, 80.0, False),
    "lag_frac": (0.0, 1.0, False),
    "gust_weight": (0.0, 1.5, False),
    "gust_decay": (5.0, 100.0, False),
}
FIXED_V2 = {**FIXED, "osc_mode": "bayes", "gust_lag": 1}
if "--space" in sys.argv and sys.argv[sys.argv.index("--space") + 1] == "v2":
    SPACE, FIXED = SPACE_V2, FIXED_V2


def run_one(args):
    params, regime, seed = args
    r = run_episode(ForecastPlanner({**params, **FIXED}), seed, regime=regime)
    return r["score"]


def evaluate(pool, params_list):
    jobs = [(p, reg, s) for p in params_list for reg, s in TRAIN_EPISODES]
    scores = pool.map(run_one, jobs, chunksize=1)
    n = len(TRAIN_EPISODES)
    return [float(np.mean(scores[i * n:(i + 1) * n])) for i in range(len(params_list))]


def to_unit(p):
    return np.array([(p[k] - lo) / (hi - lo) for k, (lo, hi, _) in SPACE.items()])


def from_unit(u):
    out = {}
    for x, (k, (lo, hi, integer)) in zip(np.clip(u, 0, 1), SPACE.items()):
        v = lo + x * (hi - lo)
        out[k] = int(round(v)) if integer else float(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gens", type=int, default=12)
    ap.add_argument("--pop", type=int, default=8)
    ap.add_argument("--sigma", type=float, default=0.25)
    ap.add_argument("--out", default=os.path.join(HERE, "params.json"))
    ap.add_argument("--space", default="v1", choices=["v1", "v2"])
    ap.add_argument("--init", default=None, help="start from these parameters (JSON) instead of the defaults")
    ap.add_argument("--log", default=os.path.join(HERE, "train_log.jsonl"))
    a = ap.parse_args()
    rng = np.random.default_rng(0)
    start = {k: DEFAULT_PARAMS[k] for k in SPACE}
    if a.init:
        start.update({k: v for k, v in json.load(open(a.init)).items() if k in SPACE})
    with Pool(os.cpu_count()) as pool:
        best_u = to_unit(start)
        best_score = evaluate(pool, [from_unit(best_u)])[0]
        sigma = a.sigma
        log = open(a.log, "w")
        print(f"gen 0  start {best_score:.3f}", flush=True)
        log.write(json.dumps({"gen": 0, "best": best_score, "params": from_unit(best_u)}) + "\n"); log.flush()
        for gen in range(1, a.gens + 1):
            t0 = time.time()
            cands = [np.clip(best_u + sigma * rng.standard_normal(len(SPACE)), 0, 1) for _ in range(a.pop)]
            scores = evaluate(pool, [from_unit(u) for u in cands])
            k = int(np.argmax(scores))
            improved = scores[k] > best_score
            if improved:
                best_u, best_score = cands[k], scores[k]
                sigma = min(sigma * 1.2, 0.4)
            else:
                sigma *= 0.7
            print(f"gen {gen}  best {best_score:.3f}  gen-max {scores[k]:.3f}  sigma {sigma:.3f}  {time.time()-t0:.0f}s", flush=True)
            log.write(json.dumps({"gen": gen, "best": best_score, "gen_scores": scores, "sigma": sigma,
                                  "params": from_unit(best_u)}) + "\n"); log.flush()
            json.dump({**from_unit(best_u), **FIXED, "train_score": best_score}, open(a.out, "w"), indent=2)


if __name__ == "__main__":
    main()
