"""A* with the TRUE future wind (an oracle, for reference): python agent/astar_oracle.py <seed> <weight> [regime]"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import SailingEnv
from forecast_agent import ForecastPlanner
from astar import astar


def run(seed, weight, regime="train", max_expansions=40000, period=5):
    env = SailingEnv(regime)
    obs, info = env.reset(seed=seed)
    helper = ForecastPlanner()
    helper.reset(info["env_info"])
    fields, maps = {}, {}

    def wind_at(t):
        if t not in fields:
            fields[t] = env.wind.field(t)
        return fields[t]

    def h(p, t):
        k = int(t // period) * period
        if k not in maps:
            maps[k] = helper._ttg_map(wind_at(k))
        return helper._value(maps[k], p)
    return astar(env.position, env.velocity, 0, wind_at, env.island_mask, h, weight=weight, max_expansions=max_expansions)


if __name__ == "__main__":
    seed, weight = int(sys.argv[1]), float(sys.argv[2])
    regime = sys.argv[3] if len(sys.argv) > 3 else "train"
    t0 = time.time()
    r = run(seed, weight, regime)
    # check the route in the real environment
    env = SailingEnv(regime); env.reset(seed=seed)
    for a in r.actions:
        _, _, term, _, _ = env.step(a)
        if term:
            break
    print(json.dumps({"seed": seed, "weight": weight, "astar": r.arrival_time, "replayed": env.arrival_time,
                      "expansions": r.expansions, "seconds": round(time.time() - t0, 1), "actions": r.actions}))
