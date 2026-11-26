"""
Observation -> network input, using only the present and the past (what a real agent sees).

Grid input (11 x 32 x 32, the map pooled by 4):
  wind now (u, v); change over the last 3 steps; change over the last 10 steps; change since step 0 (a turning front);
  island fraction; boat position (a small Gaussian blob); distance to the goal around the islands.
Vector input (44): position, velocity, step, local wind, and for each of the 9 actions the exact outcome of one step
in the present wind (displacement, crash, arrival, gain in distance to the goal).
"""
import numpy as np

from sailing_challenge import distance_to_goal
from sailing_challenge.safety import next_step

G = 32
POOL = 128 // G


def pool(a):
    return a.reshape(G, POOL, G, POOL, *a.shape[2:]).mean(axis=(1, 3))


class FeatureBuilder:
    def __init__(self, env_info):
        self.info = env_info
        self.mask = np.asarray(env_info["island_mask"], dtype=bool)
        self.D = distance_to_goal(env_info)
        Dc = np.where(np.isfinite(self.D), self.D, 200.0)
        self.mask_c = pool(self.mask.astype(np.float32))
        self.dist_c = pool(Dc.astype(np.float32)) / 128.0
        self.hist = []
        yy, xx = np.mgrid[0:G, 0:G]
        self._yy, self._xx = yy + 0.5, xx + 0.5

    def __call__(self, obs):
        f = pool(np.asarray(obs["wind_field"], dtype=np.float32)) / 10.0   # (32, 32, 2)
        self.hist.append(f)
        h = self.hist
        d3 = f - h[max(0, len(h) - 4)]
        d10 = f - h[max(0, len(h) - 11)]
        d0 = f - h[0]
        p = np.asarray(obs["position"], dtype=float)
        blob = np.exp(-((self._xx - p[0] / POOL) ** 2 + (self._yy - p[1] / POOL) ** 2) / 2.0).astype(np.float32)
        grid = np.concatenate([f, d3, d10, d0], axis=2).transpose(2, 0, 1)
        grid = np.concatenate([grid, self.mask_c[None], blob[None], self.dist_c[None]], axis=0).astype(np.float32)
        v = np.asarray(obs["velocity"], dtype=float)
        lw = np.asarray(obs["local_wind"], dtype=float)
        d_now = self._dist(p)
        vec = [p[0] / 128, p[1] / 128, v[0] / 8, v[1] / 8, obs["step"] / 100.0, lw[0] / 10, lw[1] / 10, 0.0]
        safe = np.zeros(9, dtype=bool)
        for a in range(9):
            tr = next_step(obs, a, self.info)
            dp = tr.position - p
            crash = tr.event == "crash"
            arrive = tr.event == "arrive"
            safe[a] = not crash
            gain = (d_now - (0.0 if arrive else self._dist(tr.position))) / 4.0
            vec += [dp[0] / 4, dp[1] / 4, float(crash), float(arrive) + gain]
        return grid, np.asarray(vec, dtype=np.float32), safe

    def _dist(self, p):
        d = self.D[min(int(p[1]), 127), min(int(p[0]), 127)]
        return float(d) if np.isfinite(d) else 200.0
