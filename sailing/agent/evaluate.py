"""
Test on episodes never used for training: unseen seeds of the training regime, unseen random regimes (other wind
settings and island layouts), and the 5 scenarios of the game. Compared agents:
  Trained       ForecastPlanner with the parameters of params.json
  Untrained     ForecastPlanner with its hand-set default parameters
  Persistence   the same planner, but the wind is assumed to stay as it is now (no forecast)
  Greedy        at every step, the move that lands closest to the goal (distance around the islands)
  Oracle        beam search that KNOWS the future wind (a reference, not an agent: unreachable in the real game)

    python agent/evaluate.py --out agent/results.json
"""
import argparse, json, os, sys, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from sailing_challenge import SailingEnv, run_episode, distance_to_goal, episode_score
from sailing_challenge.agents import BaseAgent
from sailing_challenge.safety import next_step
from sailing_challenge.dynamics import transition, local_wind
from forecast_agent import ForecastPlanner

TEST_EPISODES = ([("scenario", s) for s in (7, 17, 3, 30, 40)] +
                 [("train", 60000 + i) for i in range(30)] +
                 [("random", 50000 + i) for i in range(40)])
PERSIST = {"osc_weight": 0.0, "trend_weight": 0.0, "front_weight": 0.0, "lag_frac": 0.0, "horizon": 5, "beam": 16}


class Greedy(BaseAgent):
    def reset(self, env_info):
        self.info, self.D = env_info, distance_to_goal(env_info)

    def act(self, obs):
        best = None
        for a in range(9):
            tr = next_step(obs, a, self.info)
            if tr.event == "arrive":
                return a
            d = self.D[min(int(tr.position[1]), 127), min(int(tr.position[0]), 127)] if tr.event == "move" else 1e9
            if best is None or d < best[0]:
                best = (d, a)
        return best[1]


def oracle(regime, seed, B=300):
    """Beam search with the true future wind (the ForecastPlanner map, recomputed on the true field every 5 steps)."""
    env = SailingEnv(regime)
    obs, info = env.reset(seed=seed)
    helper = ForecastPlanner()
    helper.reset(info["env_info"])
    mask = env.island_mask
    fields = [env.wind.field(t) for t in range(200)]
    maps = {}

    def h(p, t):
        k = min(int(t // 5) * 5, 195)
        if k not in maps:
            maps[k] = helper._ttg_map(fields[k])
        return helper._value(maps[k], p)
    beam, best = [(np.array(env.position), np.zeros(2))], None
    for t in range(199):
        cand = []
        for p, v in beam:
            w = local_wind(fields[t], p)
            for a in range(9):
                tr = transition(p, v, w, a, mask)
                if tr.event == "arrive":
                    if best is None or t + tr.fraction < best:
                        best = t + tr.fraction
                elif tr.event == "move":
                    cand.append((t + 1 + h(tr.position, t + 1), tr.position, tr.velocity))
        if not cand or (best is not None and min(c[0] for c in cand) > best + 3):
            break
        cand.sort(key=lambda c: c[0])
        seen, beam = set(), []
        for f, p, v in cand:
            k = (round(p[0] * 2), round(p[1] * 2), round(v[0] * 2), round(v[1] * 2))
            if k not in seen:
                seen.add(k)
                beam.append((p, v))
                if len(beam) >= B:
                    break
    return best


def job(args):
    name, kind, seed, params = args
    regime = "train" if kind == "scenario" else kind
    t0 = time.time()
    if name == "Oracle":
        at = oracle(regime, seed)
        return {"agent": name, "kind": kind, "seed": seed, "arrival_time": at, "score": episode_score(at),
                "seconds": time.time() - t0}
    if name.startswith("Net:"):   # a trained network (alphastar/), e.g. "Net:alphastar/bc.pt"
        sys.path.insert(0, os.path.join(os.path.dirname(HERE), "alphastar"))
        from model import NetAgent
        agent = NetAgent(os.path.join(os.path.dirname(HERE), name[4:]))
    else:
        agent = Greedy() if name == "Greedy" else ForecastPlanner(params)
    r = run_episode(agent, seed, regime=regime)
    return {"agent": name, "kind": kind, "seed": seed, "arrival_time": r["arrival_time"], "score": r["score"],
            "crashes": r["crashes"], "seconds": time.time() - t0, "act_seconds": r["act_seconds"], "steps": r["steps"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default=os.path.join(HERE, "params.json"))
    ap.add_argument("--out", default=os.path.join(HERE, "results.json"))
    ap.add_argument("--agents", default="Trained,Untrained,Persistence,Greedy,Oracle")
    a = ap.parse_args()
    trained = {k: v for k, v in json.load(open(a.params)).items() if k != "train_score"}
    params = {"Trained": trained, "Untrained": {"horizon": 5, "beam": 16}, "Persistence": PERSIST, "Greedy": None, "Oracle": None}
    names = a.agents.split(",")
    params.update({n: None for n in names if n.startswith("Net:")})
    jobs = [(n, k, s, params[n]) for n in names for k, s in TEST_EPISODES]
    jobs.sort(key=lambda j: j[0] != "Oracle")  # the long jobs first
    with Pool(os.cpu_count()) as pool:
        res = pool.map(job, jobs, chunksize=1)
    json.dump(res, open(a.out, "w"), indent=1)
    print(f"{'agent':12s} {'scenarios':>10s} {'train-regime':>13s} {'random-regime':>14s} {'all':>7s} {'arrived':>8s}")
    for n in names:
        rows = [r for r in res if r["agent"] == n]
        m = lambda k: np.mean([r["score"] for r in rows if r["kind"] == k])
        arr = np.mean([r["arrival_time"] is not None for r in rows])
        print(f"{n:12s} {m('scenario'):10.2f} {m('train'):13.2f} {m('random'):14.2f} {np.mean([r['score'] for r in rows]):7.2f} {arr:8.0%}")


if __name__ == "__main__":
    main()
