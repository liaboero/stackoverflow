"""
Diagnostic: which part of the wind forecast costs the most? The planner is given the TRUE future of one component
of the wind (oscillation, turning front, gusts) and the present value of the others, frozen. The gain over
"everything frozen" measures what a perfect forecast of that component would be worth.

    python agent/diagnose.py
"""
import json, os, sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import SailingEnv, run_episode
from sailing_challenge.wind import rotate
from forecast_agent import ForecastPlanner

BASE = {k: v for k, v in json.load(open(os.path.join(HERE, "trained_params.json"))).items() if k != "train_score"}


class ComponentOracle(ForecastPlanner):
    """Forecast = true future for the components in `known`, present value for the others."""

    def __init__(self, wind, known, params, hybrid=False):
        super().__init__(params)
        self.wind, self.known, self.hybrid = wind, set(known), hybrid

    def _parts(self, t):
        W = self.wind
        osc = W.delta(t) if W.oscillating else 0.0
        front = 0.0
        if W.front is not None and W.front["kind"] == "turning":
            front = W.front_twist(t) - W.front_twist(0)
        shift, factor = 0.0, 1.0
        if W.gusting:
            factor, shift = W.gust_effect(t)
        return osc, front, shift, factor

    def _forecast_field(self, field, t, tau, T):
        if tau <= 0:
            return np.asarray(field, dtype=float)
        now, fut = self._parts(t), self._parts(t + int(round(tau)))
        if self.hybrid:   # the agent's own forecast for the components it does not know
            # the agent works with mathematical angles (counter-clockwise); rotate() turns the wind-from direction
            # clockwise, which turns the vector clockwise: a compass turn of -delta
            dg = self._forecast_g(t, tau)
            ds = self._forecast_s(t, tau, T)
            now = (now[0] - dg, now[1] - ds * self.lat, now[2], now[3])
        pick = lambda i, name: fut[i] if name in self.known else now[i]
        turn = pick(0, "osc") + pick(1, "front") + pick(2, "gusts")
        f = rotate(self.wind._base64, turn)
        factor = pick(3, "gusts")
        return f * np.asarray(factor, dtype=float)[..., None] if np.ndim(factor) else f * factor


HYBRID = {   # the current agent, with the true future of one component
    "agent, recomposed": [],
    "agent + true osc": ["osc"],
    "agent + true front": ["front"],
    "agent + true gusts": ["gusts"],
}
VARIANTS = {
    "current agent": None,
    "all frozen": [],
    "true oscillation": ["osc"],
    "true front": ["front"],
    "true gusts": ["gusts"],
    "all true": ["osc", "front", "gusts"],
}


def job(args):
    name, regime, seed = args
    if VARIANTS[name] is None:
        agent = ForecastPlanner(BASE)
    else:
        env = SailingEnv(regime); env.reset(seed=seed)
        agent = ComponentOracle(env.wind, VARIANTS[name], BASE, hybrid=name in HYBRID)
    return name, regime, seed, run_episode(agent, seed, regime=regime)["score"]


if __name__ == "__main__":
    if "--hybrid" in sys.argv:
        VARIANTS = {"current agent": None, **HYBRID, "all true": ["osc", "front", "gusts"]}
    eps = [("train", 90000 + i) for i in range(24)] + [("random", 91000 + i) for i in range(24)]
    with Pool(os.cpu_count()) as p:
        res = p.map(job, [(n, r, s) for n in VARIANTS for r, s in eps], chunksize=1)
    sc = {(n, r, s): x for n, r, s, x in res}
    base = "current agent" if "--hybrid" in sys.argv else "all frozen"
    ref = np.array([sc[(base, r, s)] for r, s in eps])
    print(f"{'variant':20s} {'train':>7s} {'random':>7s} {'all':>7s}   gain vs {base}")
    for n in VARIANTS:
        v = np.array([sc[(n, r, s)] for r, s in eps])
        d = v - ref
        print(f"{n:20s} {v[:24].mean():7.2f} {v[24:].mean():7.2f} {v.mean():7.2f}   {d.mean():+.2f} ± {d.std(ddof=1)/np.sqrt(len(d)):.2f}")
