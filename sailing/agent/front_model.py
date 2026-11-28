"""
Probability that the wind holds a TURNING front, learned with a direct signal: the generator knows whether each map
has one. Only the wind is needed (no sailing), so thousands of maps are cheap. One logistic regression per step
t = 0..T_MAX, on what the agent has measured so far (lateral twist of the first field, its drift since).

    python agent/front_model.py          # builds the data, fits, compares with the old threshold rule, saves front_model.json
"""
import json, math, os, sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]

T_MAX = 24


def features(agent, t):
    """Features of the agent's wind analysis at step t (after _analyse of the field of step t)."""
    c_lo, c_hi, c0 = agent.c_lo, agent.c_hi, agent.c0
    s = np.asarray(agent.s_hist)
    sg = math.copysign(1.0, c0) if c0 != 0 else 1.0
    prog = -s[-1] * sg                               # > 0: the twist moves back, as a turning front does
    slope = -(s[-1] - s[max(0, len(s) - 3)]) * sg / max(1, min(2, len(s) - 1))
    return [1.0, abs(c0) / 50, (c_lo * c_hi > 0) * 1.0, min(abs(c_lo), abs(c_hi)) / 50,
            abs(c_lo - c_hi) / 50, prog / 20, slope / 5, abs(s[-1]) / 20, abs(np.asarray(agent.g_hist)[-1]) / 30]


def episode(args):
    regime, seed = args
    from sailing_challenge import SailingEnv
    from forecast_agent import ForecastPlanner
    env = SailingEnv(regime); _, info = env.reset(seed=seed)
    f = env.wind.front
    y = int(f is not None and f["kind"] == "turning")
    a = ForecastPlanner(); a.reset(info["env_info"])
    X, old = [], []
    for t in range(T_MAX + 1):
        a._analyse(env.wind.field(t))
        X.append(features(a, t))
        a._front_presumed = False
        old.append(int(a._front_T(t) is not None))   # the old rule: a turning front is assumed or seen
    return y, X, old


def fit_logistic(X, y, l2=1e-2, iters=50):
    w = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ w))
        g = X.T @ (p - y) + l2 * w
        H = (X * (p * (1 - p))[:, None]).T @ X + l2 * np.eye(len(w))
        w -= np.linalg.solve(H, g)
    return w


if __name__ == "__main__":
    train = [("random", 600000 + i) for i in range(1500)] + [("train", 700000 + i) for i in range(1500)]
    test = [("random", 610000 + i) for i in range(500)] + [("train", 710000 + i) for i in range(500)]
    with Pool(os.cpu_count()) as p:
        tr = p.map(episode, train, chunksize=8)
        te = p.map(episode, test, chunksize=8)
    W = []
    print(" t   base-rate  logloss: model  old rule   accuracy: model  old rule")
    for t in range(T_MAX + 1):
        X = np.array([x[1][t] for x in tr]); y = np.array([x[0] for x in tr], float)
        w = fit_logistic(X, y)
        W.append(w.tolist())
        Xt = np.array([x[1][t] for x in te]); yt = np.array([x[0] for x in te], float)
        pt = np.clip(1 / (1 + np.exp(-Xt @ w)), 1e-4, 1 - 1e-4)
        old = np.array([x[2][t] for x in te], float)
        po = np.clip(old * 0.9 + (1 - old) * 0.1, 1e-4, 1 - 1e-4)   # the old rule, read as a 90% / 10% belief
        ll = lambda q: float(-(yt * np.log(q) + (1 - yt) * np.log(1 - q)).mean())
        if t % 3 == 0 or t == T_MAX:
            print(f"{t:2d}   {yt.mean():.2f}       {ll(pt):.3f}     {ll(po):.3f}       {((pt > .5) == yt).mean():.3f}         {(old == yt).mean():.3f}")
    json.dump({"T_MAX": T_MAX, "weights": W}, open(os.path.join(HERE, "front_model.json"), "w"))
