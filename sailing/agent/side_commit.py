"""
The side model at the moment the agent commits: the agent sails each labelled episode; when it is 6 rows below the
corridors, its own analysis gives the features. Is a learned model then better than the agent's own choice?

    python agent/side_commit.py
"""
import json, os, sys
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import SailingEnv
from forecast_agent import ForecastPlanner
import side_model as sm
from front_model import fit_logistic

P = json.load(open(os.path.join(HERE, "params_s1_prob.json")))


def go(r):
    env = SailingEnv(r["regime"]); obs, info = env.reset(seed=r["seed"])
    y0, y1 = sm.corridor_rows(env.island_mask)
    a = ForecastPlanner(P); a.reset(info["env_info"])
    feats, side, done, t = None, None, False, 0
    while not done:
        if feats is None and obs["position"][1] >= y0 - 6:
            a._pos = obs["position"]
            a_copy_state = (list(a.g_hist), list(a.s_hist))
            feats = sm.side_features(a, obs["wind_field"], t, 35.0) + [obs["position"][0] / 64 - 1, t / 20]
        act = a.act(obs)
        obs, _, term, trunc, _ = env.step(act)
        t += 1
        done = term or trunc
        if side is None and obs["position"][1] >= (y0 + y1) / 2:
            side = obs["position"][0] > 64.5
    return feats, bool(side)


if __name__ == "__main__":
    rows = [json.loads(l) for l in open(os.path.join(HERE, "side_labels.jsonl"))]
    with Pool(os.cpu_count()) as p:
        out = p.map(go, rows, chunksize=4)
    keep = [i for i, o in enumerate(out) if o[0] is not None]
    X = np.array([out[i][0] for i in keep]); agent_right = np.array([out[i][1] for i in keep])
    y = np.array([rows[i]["t_right"] < rows[i]["t_left"] for i in keep], float)
    m = np.array([abs(rows[i]["t_right"] - rows[i]["t_left"]) for i in keep])
    accs, accs_agent = [], []
    for fold in range(5):   # 5-fold cross-validation: every episode is predicted by a model that did not see it
        te = np.arange(len(y)) % 5 == fold
        w = fit_logistic(X[~te], y[~te], l2=1.0)
        pr = 1 / (1 + np.exp(-X[te] @ w))
        accs += list((pr > .5) == y[te]); accs_agent += list(agent_right[te] == y[te])
    accs, accs_agent = np.array(accs), np.array(accs_agent)
    order = np.concatenate([np.where(np.arange(len(y)) % 5 == f)[0] for f in range(5)])
    mm = m[order]
    w_all = fit_logistic(X, y, l2=1.0)
    json.dump({"weights": w_all.tolist(), "rows_below": 6}, open(os.path.join(HERE, "side_model.json"), "w"))
    print(f"episodes {len(y)}: side model at the commit point {accs.mean():.1%} (weighted {np.average(accs, weights=mm):.1%})"
          f"   agent's own choice {accs_agent.mean():.1%} (weighted {np.average(accs_agent, weights=mm):.1%})")
