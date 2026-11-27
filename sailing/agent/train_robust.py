"""
Training that resists over-fitting: every generation sails FRESH episodes (3/4 random regimes, 1/4 training regime),
and the best candidate replaces the incumbent only if it beats it again on a second batch of fresh episodes (two
stages: the first batch selects, the second one confirms without the selection bias). The seeds never meet the
training sets of train.py, the development set, the test set or the validation set.

    python agent/train_robust.py --init agent/best_params.json --gens 10 --pop 6 --out agent/params_robust.json
"""
import argparse, json, os, sys, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import run_episode
from forecast_agent import ForecastPlanner
from train import SPACE_V2 as SPACE, FIXED_V2 as FIXED, to_unit as _tu, from_unit as _fu
import train

train.SPACE = SPACE   # to_unit / from_unit read the module's SPACE


def batch(gen, stage, n=48):
    k = 3 * n // 4
    base = 400000 + 1000 * gen + 500 * stage
    return [("random", base + i) for i in range(k)] + [("train", base + 100000 + i) for i in range(n - k)]


def job(args):
    params, regime, seed = args
    return run_episode(ForecastPlanner({**params, **FIXED}), seed, regime=regime)["score"]


def evaluate(pool, params_list, episodes):
    res = pool.map(job, [(p, r, s) for p in params_list for r, s in episodes], chunksize=1)
    n = len(episodes)
    return [np.array(res[i * n:(i + 1) * n]) for i in range(len(params_list))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--init", required=True)
    ap.add_argument("--gens", type=int, default=10)
    ap.add_argument("--pop", type=int, default=6)
    ap.add_argument("--sigma", type=float, default=0.15)
    ap.add_argument("--out", default=os.path.join(HERE, "params_robust.json"))
    ap.add_argument("--log", default=os.path.join(HERE, "train_robust_log.jsonl"))
    a = ap.parse_args()
    rng = np.random.default_rng(1)
    init = json.load(open(a.init))
    inc = {k: init[k] for k in SPACE}
    sigma = a.sigma
    log = open(a.log, "w")
    with Pool(os.cpu_count()) as pool:
        for gen in range(a.gens):
            t0 = time.time()
            u = _tu(inc)
            cands = [_fu(np.clip(u + sigma * rng.standard_normal(len(SPACE)), 0, 1)) for _ in range(a.pop)]
            eps = batch(gen, 0)
            sc = evaluate(pool, [inc] + cands, eps)
            gains = [float((s - sc[0]).mean()) for s in sc[1:]]
            k = int(np.argmax(gains))
            accepted, confirm = False, None
            if gains[k] > 0:
                eps2 = batch(gen, 1)
                s_inc, s_c = evaluate(pool, [inc, cands[k]], eps2)
                d = s_c - s_inc
                confirm = (float(d.mean()), float(d.std(ddof=1) / np.sqrt(len(d))))
                accepted = confirm[0] > 0
            if accepted:
                inc = cands[k]
                sigma = min(sigma * 1.2, 0.4)
            else:
                sigma = max(sigma * 0.8, 0.03)
            msg = (f"gen {gen + 1}  selection gain {gains[k]:+.2f}  confirmation "
                   + (f"{confirm[0]:+.2f} ± {confirm[1]:.2f}" if confirm else "-")
                   + f"  {'ACCEPTED' if accepted else 'rejected'}  sigma {sigma:.3f}  {time.time() - t0:.0f}s")
            print(msg, flush=True)
            log.write(json.dumps({"gen": gen + 1, "gains": gains, "confirm": confirm, "accepted": accepted,
                                  "params": inc, "sigma": sigma}) + "\n"); log.flush()
            json.dump({**inc, **FIXED}, open(a.out, "w"), indent=2)


if __name__ == "__main__":
    main()
