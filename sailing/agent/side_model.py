"""
Which side of the central island is faster? A direct signal for the side choice.

Labels: the planner with the TRUE future wind sails each episode twice, once with the left corridor closed and once
with the right one closed; the faster side is the label, the time difference its weight.
Features: what the agent knows at step 0 (its own expected times by each side under its forecasts, the front
probability, the mean wind, the island layout).

    python agent/side_model.py label  N   # labels N episodes per regime (seeds 800000+ / 900000+)
    python agent/side_model.py fit        # fits the side model and reports its accuracy on held-out episodes
"""
import json, math, os, sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import SailingEnv, run_episode
from forecast_agent import ForecastPlanner
from diagnose import ComponentOracle

BEST = {k: v for k, v in json.load(open(os.path.join(HERE, "best_params.json"))).items() if k != "train_score"}
LABELS = os.path.join(HERE, "side_labels.jsonl")


def corridor_rows(mask):
    """Rows where the central island stands (its columns around x = 64)."""
    rows = np.where(mask[:, 58:71].any(axis=1))[0]
    if not len(rows):
        return 0, -1
    end = rows[0]
    for r in rows[1:]:          # the first run of rows only (not the small island just south of the goal)
        if r != end + 1:
            break
        end = r
    return int(rows[0]), int(end)


def corridor(mask, side):
    """Cells of one corridor: beside the central island, on the left (x < 64.5) or the right side."""
    y0, y1 = corridor_rows(mask)
    c = np.zeros_like(mask)
    if y1 >= y0:
        if side == "left":
            c[y0:y1 + 1, :64] = True
        else:
            c[y0:y1 + 1, 65:] = True
    return c & ~mask


class ClosedSide(ComponentOracle):
    """The planner with the true future wind, which may not use one corridor."""

    def __init__(self, wind, closed, params):
        super().__init__(wind, ["osc", "front", "gusts"], params)
        self.closed = closed

    def reset(self, env_info):
        super().reset(env_info)
        self.free = self.free & ~corridor(self.mask, self.closed)


def label(args):
    regime, seed = args
    out = {"regime": regime, "seed": seed}
    for closed, key in (("right", "t_left"), ("left", "t_right")):
        env = SailingEnv(regime); env.reset(seed=seed)
        r = run_episode(ClosedSide(env.wind, closed, BEST), seed, regime=regime)
        out[key] = r["arrival_time"] if r["arrival_time"] is not None else 500.0
    return out


def from_map(agent, field, src):
    """Time from the boat's position to every cell (forward Dijkstra, the same edge times as _ttg_map)."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra
    from sailing_challenge.physics import target_velocities
    from forecast_agent import MOVES, N
    tv = target_velocities(field)[:8]
    speed = np.maximum(np.hypot(tv[..., 0], tv[..., 1]), agent.p["speed_floor"])
    idx = np.arange(N * N).reshape(N, N)
    rows, cols, data = [], [], []
    ys, xs = np.mgrid[0:N, 0:N]
    for a, (dx, dy) in enumerate(MOVES):
        u, v = xs + dx, ys + dy
        ok = (u >= 0) & (u < N) & (v >= 0) & (v < N) & agent.free
        ok[ok] &= agent.free[v[ok], u[ok]]
        rows.append(idx[ys[ok], xs[ok]]); cols.append(idx[v[ok], u[ok]])
        data.append(math.hypot(dx, dy) / speed[a][ok])
    G = csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(N * N, N * N))
    return dijkstra(G, directed=True, indices=int(src[1]) * N + int(src[0])).reshape(N, N)


def side_features(agent, field, t, lag):
    """Features of the side choice, from the agent's own analysis at step t."""
    p_turn, T_turn = agent._front_prob_T(t)
    left, right = corridor(agent.mask, "left"), corridor(agent.mask, "right")

    y0, y1 = corridor_rows(agent.mask)
    mid = np.zeros_like(agent.mask)
    mid[(y0 + y1) // 2, :] = True          # the middle row of the corridors
    pos = agent._pos

    def via(F, Tm, c):
        """Fastest route through one corridor: from the boat to a cell of its middle row, then to the goal."""
        tot = F + Tm
        v = tot[c & mid & np.isfinite(tot)]
        return float(v.min()) if len(v) else 300.0
    d = []
    for T in (T_turn, None):
        fut = agent._forecast_field(field, t, lag, T)
        Tm, F = agent._ttg_map(fut), from_map(agent, agent._forecast_field(field, t, lag / 3, T), pos)
        d.append((via(F, Tm, right) - via(F, Tm, left)) / 10)
    Tnow, Fnow = agent._ttg_map(field), from_map(agent, field, pos)
    d_now = (via(Fnow, Tnow, right) - via(Fnow, Tnow, left)) / 10
    ang = math.radians(float(np.degrees(np.arctan2(field[..., 1].mean(), field[..., 0].mean()))))
    lay = agent.mask.sum() / 1000
    return [1.0, d[0] * p_turn + d[1] * (1 - p_turn), d[0], d[1], d_now, p_turn, (2 * p_turn - 1) * math.copysign(1, agent.c0),
            math.cos(ang), math.sin(ang), lay]


def features_episode(args):
    regime, seed = args
    env = SailingEnv(regime); obs, info = env.reset(seed=seed)
    a = ForecastPlanner({**BEST, "front_mode": "prob"}); a.reset(info["env_info"])
    a._analyse(obs["wind_field"])
    a._pos = obs["position"]
    return side_features(a, obs["wind_field"], 0, 35.0)


if __name__ == "__main__":
    if sys.argv[1] == "label":
        n = int(sys.argv[2])
        eps = [("random", 800000 + i) for i in range(n)] + [("train", 900000 + i) for i in range(n)]
        with Pool(os.cpu_count()) as p, open(LABELS, "a") as f:
            for r in p.imap_unordered(label, eps, chunksize=1):
                f.write(json.dumps(r) + "\n"); f.flush()
    elif sys.argv[1] == "fit":
        from front_model import fit_logistic
        rows = [json.loads(l) for l in open(LABELS)]
        with Pool(os.cpu_count()) as p:
            X = np.array(p.map(features_episode, [(r["regime"], r["seed"]) for r in rows]))
        y = np.array([r["t_right"] < r["t_left"] for r in rows], float)
        margin = np.array([abs(r["t_right"] - r["t_left"]) for r in rows])
        wts = np.minimum(margin, 20.0) + 0.5   # a side that matters more counts more
        n = len(rows); k = int(n * 0.75)
        idx = np.random.default_rng(0).permutation(n); tr, te = idx[:k], idx[k:]
        # weighted logistic regression (duplicate by weight through the gradient)
        w = np.zeros(X.shape[1])
        for _ in range(60):
            pr = 1 / (1 + np.exp(-X[tr] @ w))
            g = X[tr].T @ ((pr - y[tr]) * wts[tr]) + 1e-2 * w
            H = (X[tr] * (pr * (1 - pr) * wts[tr])[:, None]).T @ X[tr] + 1e-2 * np.eye(len(w))
            w -= np.linalg.solve(H, g)
        pt = 1 / (1 + np.exp(-X[te] @ w))
        acc = ((pt > .5) == y[te])
        own = ((X[te][:, 1] < 0) == y[te])   # the agent's own expected-time difference alone
        print(f"labelled episodes {n}  right side better in {y.mean():.0%}  median margin {np.median(margin):.1f} steps")
        print(f"held-out accuracy: side model {acc.mean():.1%}  (weighted by the margin: {np.average(acc, weights=margin[te]):.1%})"
              f"   agent's own map alone {own.mean():.1%} (weighted {np.average(own, weights=margin[te]):.1%})")
        json.dump({"weights": w.tolist()}, open(os.path.join(HERE, "side_model.json"), "w"))
