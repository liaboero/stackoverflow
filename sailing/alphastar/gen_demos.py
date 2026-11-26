"""
Expert demonstrations: the Oracle (beam search that knows the future wind) sails each episode; every state it goes
through is stored with the features the student will see (present and past wind only) and the Oracle's action.

    python alphastar/gen_demos.py --n 240 --out /path/demos.npz
"""
import argparse, os, sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [ROOT, HERE, os.path.join(ROOT, "agent")]
from sailing_challenge import SailingEnv
from sailing_challenge.dynamics import transition, local_wind
from forecast_agent import ForecastPlanner
from features import FeatureBuilder


def oracle_actions(env, fields, start_t=0, pos=None, vel=None, B=200):
    helper = ForecastPlanner()
    helper.reset({"island_mask": env.island_mask, "goal": (64.5, 127.5), "goal_radius": 2.0})
    mask, maps = env.island_mask, {}

    def h(p, t):
        k = min(int(t // 5) * 5, len(fields) - 1)
        if k not in maps:
            maps[k] = helper._ttg_map(fields[k])
        return helper._value(maps[k], p)
    p0 = np.array(env.position if pos is None else pos, dtype=float)
    v0 = np.zeros(2) if vel is None else np.asarray(vel, dtype=float)
    beam, best = [(p0, v0, [])], None
    for t in range(start_t, len(fields) - 1):
        cand = []
        for p, v, acts in beam:
            w = local_wind(fields[t], p)
            for a in range(9):
                tr = transition(p, v, w, a, mask)
                if tr.event == "arrive":
                    if best is None or t + tr.fraction < best[0]:
                        best = (t + tr.fraction, acts + [a])
                elif tr.event == "move":
                    cand.append((t + 1 + h(tr.position, t + 1), tr.position, tr.velocity, acts + [a]))
        if not cand or (best is not None and min(c[0] for c in cand) > best[0] + 3):
            break
        cand.sort(key=lambda c: c[0])
        seen, beam = set(), []
        for f, p, v, acts in cand:
            k = (round(p[0] * 2), round(p[1] * 2), round(v[0] * 2), round(v[1] * 2))
            if k not in seen:
                seen.add(k)
                beam.append((p, v, acts))
                if len(beam) >= B:
                    break
    return best


def episode(args):
    regime, seed = args
    env = SailingEnv(regime)
    obs, info = env.reset(seed=seed)
    fields = [env.wind.field(t) for t in range(220)]
    best = oracle_actions(env, fields)
    if best is None:
        return None
    fb = FeatureBuilder(info["env_info"])
    grids, vecs, safes, acts = [], [], [], []
    for a in best[1]:
        g, v, s = fb(obs)
        grids.append(g.astype(np.float16)); vecs.append(v); safes.append(s); acts.append(a)
        obs, _, term, trunc, _ = env.step(a)
        if term or trunc:
            break
    return np.stack(grids), np.stack(vecs), np.stack(safes), np.array(acts, dtype=np.int64), best[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=240)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    a = ap.parse_args()
    # seeds disjoint from training (1000+, 10000+) and test (60000+, 50000+, scenarios) of the planner
    eps = [("train" if i % 2 == 0 else "random", 200000 + i) for i in range(a.n)]
    with Pool(a.procs) as pool:
        res = [r for r in pool.imap(episode, eps, chunksize=1) if r is not None]
    G, V, S, A, T = zip(*res)
    ep_id = np.concatenate([np.full(len(x), i) for i, x in enumerate(A)])
    np.savez(a.out, grid=np.concatenate(G), vec=np.concatenate(V), safe=np.concatenate(S), act=np.concatenate(A),
             ep=ep_id, times=np.array(T))
    print("episodes", len(res), "samples", len(ep_id), "mean oracle time", float(np.mean(T)))


if __name__ == "__main__":
    main()
