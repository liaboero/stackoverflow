"""
ForecastPlanner: an agent that only sees the wind of the present and the past (never the future), forecasts how the
wind field will turn, and plans its route on that forecast.

Every step:
  1. Wind analysis. The direction of every cell is compared with the field of step 0; the difference is fitted by
       delta(x, y, t) = g(t) + s(t) * tanh((x - 64.5) / 20)
     g = global rotation (the oscillation), s = lateral twist (a front that turns, opposite on the two sides).
  2. Forecast. g is extrapolated with a sinusoid fitted on its history (period searched in a grid), s with the shape of
     a turning front, s(t) = c0 * (cos(pi * min(t / T, 1)) - 1), whose amplitude c0 is read on the field of step 0
     and whose duration T is estimated from the history (a prior until the turn is visible).
  3. Planning. A time-to-goal map (Dijkstra on the grid, edge time = length / best speed of the boat in that direction)
     is computed on the forecast field at the time the boat is expected to be half-way; a short beam search with the
     exact boat model, on the forecast fields of the next steps, ends on that map.
The weights of the forecasts, the planning horizon and the priors are the parameters learned by train.py.
"""
import math

import numpy as np

from sailing_challenge.agents import BaseAgent
from sailing_challenge.dynamics import transition
from sailing_challenge.physics import target_velocities

try:
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra as _sp_dijkstra
except ImportError:  # e.g. in the browser without scipy: a pure Python Dijkstra
    csr_matrix = _sp_dijkstra = None

N = 128
MOVES = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]

DEFAULT_PARAMS = {
    "osc_weight": 1.0,       # 0: the global rotation stays as it is now, 1: sinusoid forecast
    "osc_min_steps": 10,     # history needed before trusting the sinusoid fit
    "trend_weight": 0.5,     # before that: weight of the linear trend of g
    "front_weight": 1.0,     # weight of the turning-front forecast
    "front_threshold": 25.0, # |c0| (degrees) above which a lateral twist of the same sign at the bottom and the top
                             # of the map is taken for a turning front
    "front_prior_T": 25.0,   # duration of the turn assumed while it is not visible yet
    "lag_frac": 0.5,         # the map is computed on the forecast at t + lag_frac * (expected remaining time)
    "lag_max": 60.0,
    "horizon": 6,            # depth of the beam search
    "beam": 24,              # width of the beam search
    "speed_floor": 0.05,     # minimal speed in the map (cells per step), avoids infinite times
}


def _wrap(a):
    return (a + 180.0) % 360.0 - 180.0


class ForecastPlanner(BaseAgent):
    def __init__(self, params=None):
        self.p = dict(DEFAULT_PARAMS)
        if params:
            self.p.update(params)

    # ------------------------------------------------------------------ episode
    def reset(self, env_info):
        self.mask = np.asarray(env_info["island_mask"], dtype=bool)
        self.goal = tuple(env_info["goal"])
        self.goal_radius = float(env_info["goal_radius"])
        self.free = ~self.mask
        xs = np.arange(N) + 0.5
        self.lat = np.tanh((xs - 64.5) / 20.0)[None, :].repeat(N, 0)          # lateral basis, [y, x]
        w = self.free.astype(float)
        self._w = w / w.sum()
        self.ang0 = None
        self.g_hist, self.s_hist = [], []
        self.c0 = 0.0
        self.turning_prob = 0.0
        self._map_cache = None

    # ------------------------------------------------------------------ wind analysis
    @staticmethod
    def _angle(field):
        return np.degrees(np.arctan2(field[..., 1], field[..., 0]))

    def _analyse(self, field):
        ang = self._angle(field)
        if self.ang0 is None:
            self.ang0 = ang
            # amplitude of a lateral twist in the field of step 0, at the bottom and at the top of the map
            lo, hi = slice(0, 40), slice(85, 128)
            c_lo = self._lateral_coef(ang[lo], self.lat[lo], self.free[lo])
            c_hi = self._lateral_coef(ang[hi], self.lat[hi], self.free[hi])
            same = c_lo * c_hi > 0
            self.c0 = 0.5 * (c_lo + c_hi)
            self.turning_prior = bool(same and min(abs(c_lo), abs(c_hi)) > self.p["front_threshold"])
        d = _wrap(ang - self.ang0)
        # least squares on [1, lat], weighted by the free cells
        w = self._w
        m1, ml, mll = w.sum(), (w * self.lat).sum(), (w * self.lat ** 2).sum()
        b1, bl = (w * d).sum(), (w * d * self.lat).sum()
        det = m1 * mll - ml * ml
        g = (b1 * mll - bl * ml) / det
        s = (m1 * bl - ml * b1) / det
        self.g_hist.append(g)
        self.s_hist.append(s)

    @staticmethod
    def _lateral_coef(ang, lat, free):
        a = ang[free]
        l = lat[free]
        # circular mean of the angle, then regression of the deviation on the lateral basis
        m = math.degrees(math.atan2(np.sin(np.radians(a)).mean(), np.cos(np.radians(a)).mean()))
        dev = _wrap(a - m)
        l = l - l.mean()
        return float((dev * l).sum() / max((l * l).sum(), 1e-9))

    # ------------------------------------------------------------------ forecasts
    def _forecast_g(self, t, tau):
        """Predicted change of the global rotation between now (t) and t + tau."""
        g = np.asarray(self.g_hist)
        p = self.p
        if tau <= 0 or len(g) < 3:
            return 0.0
        if len(g) < p["osc_min_steps"]:
            k = min(len(g) - 1, 4)
            rate = (g[-1] - g[-1 - k]) / k
            return p["trend_weight"] * rate * min(tau, 15.0)
        ts = np.arange(len(g), dtype=float)
        best = None
        for P in np.arange(40.0, 101.0, 4.0):
            X = np.stack([np.sin(2 * np.pi * ts / P), np.cos(2 * np.pi * ts / P), np.ones_like(ts)], 1)
            coef, res, *_ = np.linalg.lstsq(X, g, rcond=None)
            r = float(((X @ coef - g) ** 2).sum())
            amp = math.hypot(coef[0], coef[1])
            r += 0.02 * len(g) * max(0.0, amp - 50.0) ** 2   # prior: amplitudes stay below ~45 degrees
            if best is None or r < best[0]:
                best = (r, P, coef)
        _, P, c = best
        f = lambda u: c[0] * math.sin(2 * math.pi * u / P) + c[1] * math.cos(2 * math.pi * u / P) + c[2]
        return p["osc_weight"] * (f(t + tau) - f(t))

    def _front_T(self, t):
        """Duration of the turn of the front (estimated), or None when there seems to be no turning front."""
        s = np.asarray(self.s_hist)
        c0 = self.c0
        if abs(c0) < 1e-6:
            return None
        ratio = 1.0 + s[-1] / c0 if len(s) else 1.0     # = cos(pi t / T) for a turning front
        moved = abs(s[-1]) > max(6.0, 0.15 * abs(c0)) and s[-1] * c0 < 0
        if moved and t > 0:
            r = min(max(ratio, -1.0), 1.0)
            if r <= -0.999:
                return float(t)  # the turn is over
            return float(min(max(math.pi * t / max(math.acos(r), 1e-3), t), 60.0))
        if self.turning_prior and t < self.p["front_prior_T"]:
            return float(self.p["front_prior_T"])
        return None

    def _forecast_s(self, t, tau, T):
        if T is None or tau <= 0:
            return 0.0
        c0 = self.c0
        now = c0 * (math.cos(math.pi * min(t / T, 1.0)) - 1.0)
        fut = c0 * (math.cos(math.pi * min((t + tau) / T, 1.0)) - 1.0)
        return self.p["front_weight"] * (fut - now)

    def _forecast_field(self, field, t, tau, T):
        dg = self._forecast_g(t, tau)
        ds = self._forecast_s(t, tau, T)
        if abs(dg) < 1e-9 and abs(ds) < 1e-9:
            return np.asarray(field, dtype=float)
        rot = np.radians(dg + ds * self.lat)
        c, s = np.cos(rot), np.sin(rot)
        f = np.asarray(field, dtype=float)
        out = np.empty_like(f)
        out[..., 0] = c * f[..., 0] - s * f[..., 1]
        out[..., 1] = s * f[..., 0] + c * f[..., 1]
        return out

    # ------------------------------------------------------------------ time-to-goal map
    def _ttg_map(self, field):
        tv = target_velocities(field)[:8]                          # (8, y, x, 2)
        speed = np.maximum(np.hypot(tv[..., 0], tv[..., 1]), self.p["speed_floor"])
        gx, gy = int(self.goal[0]), int(self.goal[1])
        if _sp_dijkstra is not None:
            idx = np.arange(N * N).reshape(N, N)
            rows, cols, data = [], [], []
            for a, (dx, dy) in enumerate(MOVES):
                ys, xs = np.mgrid[0:N, 0:N]
                u, v = xs + dx, ys + dy
                ok = (u >= 0) & (u < N) & (v >= 0) & (v < N)
                ok &= self.free
                ok[ok] &= self.free[v[ok], u[ok]]
                # reversed graph: edge (target -> source) so that one Dijkstra from the goal gives every time-to-goal
                rows.append(idx[v[ok], u[ok]])
                cols.append(idx[ys[ok], xs[ok]])
                data.append(math.hypot(dx, dy) / speed[a][ok])
            G = csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(N * N, N * N))
            T = _sp_dijkstra(G, directed=True, indices=gy * N + gx).reshape(N, N)
        else:
            import heapq
            T = np.full((N, N), np.inf)
            T[gy, gx] = 0.0
            heap = [(0.0, gx, gy)]
            sp = speed.tolist()
            free = self.free.tolist()
            while heap:
                d, x, y = heapq.heappop(heap)
                if d > T[y, x]:
                    continue
                for a, (dx, dy) in enumerate(MOVES):
                    u, v = x - dx, y - dy
                    if 0 <= u < N and 0 <= v < N and free[v][u]:
                        c = d + math.hypot(dx, dy) / sp[a][v][u]
                        if c < T[v, u]:
                            T[v, u] = c
                            heapq.heappush(heap, (c, u, v))
        # cells within the goal circle cost nothing more
        return np.maximum(T - self.goal_radius / 3.0, 0.0)

    @staticmethod
    def _value(T, p):
        """Bilinear interpolation of the map at a continuous position (cell centres at +0.5)."""
        x = min(max(p[0] - 0.5, 0.0), N - 1.001)
        y = min(max(p[1] - 0.5, 0.0), N - 1.001)
        i, j = int(x), int(y)
        fx, fy = x - i, y - j
        vals = T[j:j + 2, i:i + 2]
        if not np.isfinite(vals).all():
            return float(T[min(int(p[1]), N - 1), min(int(p[0]), N - 1)])
        return float((vals[0, 0] * (1 - fx) + vals[0, 1] * fx) * (1 - fy) + (vals[1, 0] * (1 - fx) + vals[1, 1] * fx) * fy)

    # ------------------------------------------------------------------ act
    def act(self, obs):
        p = self.p
        field = obs["wind_field"]
        t = int(obs["step"])
        self._analyse(field)
        T_front = self._front_T(t)
        pos = np.asarray(obs["position"], dtype=float)
        vel = np.asarray(obs["velocity"], dtype=float)

        # map on the forecast at the expected half-way time (two passes: the remaining time depends on the map)
        T0 = self._ttg_map(field) if self._map_cache is None else self._map_cache
        remaining = self._value(T0, pos)
        if not np.isfinite(remaining):
            remaining = 60.0
        lag = min(p["lag_frac"] * remaining, p["lag_max"])
        Tmap = self._ttg_map(self._forecast_field(field, t, lag, T_front))
        self._map_cache = Tmap

        # beam search on the forecast fields of the next steps
        H, B = int(p["horizon"]), int(p["beam"])
        fields = [self._forecast_field(field, t, k, T_front) for k in range(H)]
        beam = [(pos, vel, None)]
        best_arrival = None
        scored = []
        for k in range(H):
            cand = []
            for q, v, first in beam:
                w = fields[k][int(q[1]), int(q[0])]
                for a in range(9):
                    tr = transition(q, v, w, a, self.mask)
                    f0 = a if first is None else first
                    if tr.event == "arrive":
                        at = k + tr.fraction
                        if best_arrival is None or at < best_arrival[0]:
                            best_arrival = (at, f0)
                    elif tr.event == "move":
                        cand.append((k + 1 + self._value(Tmap, tr.position), tr.position, tr.velocity, f0))
            if best_arrival is not None and (not cand or best_arrival[0] <= min(c[0] for c in cand)):
                return int(best_arrival[1])
            if not cand:
                break
            cand.sort(key=lambda c: c[0])
            seen, beam = set(), []
            for f, q, v, f0 in cand:
                key = (round(q[0] * 4), round(q[1] * 4), round(v[0] * 2), round(v[1] * 2))
                if key in seen:
                    continue
                seen.add(key)
                beam.append((q, v, f0))
                if len(beam) >= B:
                    break
            scored = cand
        if best_arrival is not None:
            return int(best_arrival[1])
        if scored:
            return int(scored[0][3])
        return 8
