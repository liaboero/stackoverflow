"""
False routes on the training episodes: the current model is compared with the same planner given the TRUE future
wind ("all true"). A large gap is a forecast error, not a planning limit. Each such episode is then diagnosed by
giving the agent the true future of one component, and by switching off one of its own forecasts.

    python agent/false_routes.py
"""
import json, os, sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import SailingEnv, run_episode
from forecast_agent import ForecastPlanner
from diagnose import ComponentOracle

os.environ["SAILING_TRAIN_EPISODES"] = "64"
from train import TRAIN_EPISODES

BEST = {k: v for k, v in json.load(open(os.path.join(HERE, "best_params.json"))).items() if k != "train_score"}
VARIANTS = {
    "model": ("own", {}),
    "all true": ("oracle", ["osc", "front", "gusts"]),
    "model + true osc": ("hybrid", ["osc"]),
    "model + true front": ("hybrid", ["front"]),
    "model + true gusts": ("hybrid", ["gusts"]),
    "model, no gust forecast": ("own", {"gust_weight": 0.0}),
    "model, no front forecast": ("own", {"front_weight": 0.0}),
    "model, no osc forecast": ("own", {"osc_weight": 0.0}),
}


def job(args):
    name, regime, seed = args
    kind, arg = VARIANTS[name]
    if kind == "own":
        agent = ForecastPlanner({**BEST, **arg})
    else:
        env = SailingEnv(regime); env.reset(seed=seed)
        agent = ComponentOracle(env.wind, arg, BEST, hybrid=kind == "hybrid")
    r = run_episode(agent, seed, regime=regime, record=True)
    return name, regime, seed, r["arrival_time"], r["trajectory"]


def describe(regime, seed):
    env = SailingEnv(regime); env.reset(seed=seed)
    W = env.wind
    f = W.front
    front = "none" if f is None else f"{f['kind']} ({'right' if f['side'] > 0 else 'left'} side lures)"
    return {"islands": env.regime.get("islands"), "front": front, "osc_amp": round(W.amplitude_deg or 0, 1),
            "osc_period": round(W.period or 0, 1), "gusts": len(W.gusts)}


def side(traj, y=40.0):
    """x when the trajectory first crosses height y (which side of the central island it takes)."""
    for p in traj:
        if p[1] >= y:
            return p[0]
    return traj[-1][0]


if __name__ == "__main__":
    eps = TRAIN_EPISODES
    with Pool(os.cpu_count()) as p:
        res = p.map(job, [(n, r, s) for n in VARIANTS for r, s in eps], chunksize=1)
    R = {(n, r, s): (t, tr) for n, r, s, t, tr in res}
    rows = []
    for r, s in eps:
        t = {n: R[(n, r, s)][0] or 500.0 for n in VARIANTS}
        rows.append({"regime": r, "seed": s, **{k: round(v, 1) for k, v in t.items()},
                     "gap": round(t["model"] - t["all true"], 1),
                     "side model": round(side(R[("model", r, s)][1]), 0),
                     "side all true": round(side(R[("all true", r, s)][1]), 0), **describe(r, s)})
    json.dump(rows, open(os.path.join(HERE, "false_routes.json"), "w"), indent=1)
    gaps = np.array([x["gap"] for x in rows])
    print(f"episodes {len(rows)}  mean gap {gaps.mean():.2f} steps  median {np.median(gaps):.2f}  "
          f"gap > 5 steps: {(gaps > 5).sum()}  share of the total gap from those: {gaps[gaps > 5].sum() / gaps[gaps > 0].sum():.0%}")
    for x in sorted(rows, key=lambda x: -x["gap"])[:14]:
        print(json.dumps(x))
