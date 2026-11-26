"""
Wind. A single generator; a regime is this generator with its settings, an episode is a seed in a regime.

The wind only exists on the grid (one vector per cell) and the boat uses the wind of its cell.
It is exogenous: it does not depend on the boat, so the same episode can be replayed for several agents.

Mechanisms (added one at a time, each measured with the benchmark):
  1. base field. Per episode: the mean direction the wind comes from (mean_from_deg, drawn uniformly, or with a
     triangular law peaked in the middle of the range if mean_from_shape is "triangular") and a spatial variability S
     (spatial_range_deg). A large-scale pattern, a cubic B-spline surface with n x n random control points (like the
     "quadrants" of v2), is scaled so that the wind direction varies by exactly S degrees across the map; a smaller-scale
     smooth noise is added. Constant speed. Two neighbouring cells never differ by more than max_neighbour_deg:
     the pattern is drawn again (from the same random stream) until it fits, and S is reduced as a last resort.
  2. oscillation. The whole field turns by delta(t) = A sin(2 pi t / P + phi) + noise(t) (degrees, positive = the
     direction the wind comes from turns clockwise). P, A and phi are drawn per episode: they are hidden from the agent.
     The noise is an AR(1) process.
  3. front (optional): in some episodes, a wind front lies across the map. Below the front, the wind on one side of
     the central island is turned away from the north (it carries the boat), and on the other side towards the north
     (head wind); above the front, the other way round. The side that carries the boat at the start is the side the
     mean wind blows to. Either the front stays still (the side that looks good at the start turns bad higher up),
     or the wind of the whole lure side turns against the boat during the race (and the other side the opposite way).
  4. gusts and lulls. Smooth spots (Gaussian, radius of tens of cells) that drift with the
     wind at constant velocity, fade in and out, and change the wind speed (strength > 0: gust, < 0: lull);
     each one can also turn the wind (shift_deg, option b). Their number, size, strength and motion follow the regime.
     The strength is drawn uniformly ("strength"), or as a clear gust or lull ("strength_magnitude", see gust_strength).
With static_wind=True, the wind is the base field and never changes (for the notebooks).

Randomness: each part of the scenario has its own random stream derived from the seed, and only uniform draws
(Generator.random) are used, so that a seed gives the same scenario whatever the NumPy version.
"""
import math

import numpy as np

from .constants import GRID_SIZE

# One random stream per part of the scenario. Never renumber: it would change every scenario.
STREAM_BASE_DIRECTION = 1
STREAM_BASE_LARGE = 2
STREAM_BASE_SMALL = 3
STREAM_BASE_SPATIAL = 4
STREAM_OSCILLATION = 5
STREAM_OSCILLATION_NOISE = 6
STREAM_GUSTS = 7
# stream 8 was used by a drifting pattern, explored on 20/11 and removed: never reuse it
STREAM_FRONT = 9

MAX_PATTERN_DRAWS = 30


def rng_stream(seed, stream):
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([int(seed), int(stream)])))


def uniform(rng, low, high, size=None):
    return low + (high - low) * rng.random(size)


def triangular(u, low, high, mode=None):
    """Inverse CDF of the triangular law on [low, high] (peak at `mode`, the middle by default) at u in [0, 1)."""
    mode = 0.5 * (low + high) if mode is None else float(mode)
    if u < (mode - low) / (high - low):
        return low + math.sqrt(u * (high - low) * (mode - low))
    return high - math.sqrt((1.0 - u) * (high - low) * (high - mode))


def standard_normal(rng, size):
    """Box-Muller transform of uniform draws (stable across NumPy versions, unlike Generator.normal)."""
    u1 = 1.0 - rng.random(size)  # in (0, 1]
    u2 = rng.random(size)
    return np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * np.pi * u2)


def bspline_basis(n_control, size=GRID_SIZE):
    """
    Matrix B of shape (size, n_control + 2): row i holds the weights of the control points at cell centre i + 0.5.
    Uniform cubic B-spline with n_control control points spread over [0, size], plus one outside each end.
    """
    h = size / (n_control - 1)
    knots = (np.arange(n_control + 2) - 1) * h
    u = np.abs((np.arange(size)[:, None] + 0.5 - knots[None, :]) / h)
    return np.where(u < 1.0, (4.0 - 6.0 * u ** 2 + 3.0 * u ** 3) / 6.0,
                    np.where(u < 2.0, (2.0 - u) ** 3 / 6.0, 0.0))


def smooth_field(rng, n_control, std, size=GRID_SIZE):
    """Smooth random surface (size x size): cubic B-spline with independent N(0, std^2) control points."""
    coeffs = std * standard_normal(rng, (n_control + 2, n_control + 2))
    basis = bspline_basis(n_control, size)
    return basis @ coeffs @ basis.T


def unit_range_pattern(rng, n_control, size=GRID_SIZE):
    """Smooth random surface rescaled to the range [-0.5, 0.5] exactly."""
    g = smooth_field(rng, n_control, 1.0, size)
    lo, hi = g.min(), g.max()
    return (g - 0.5 * (lo + hi)) / (hi - lo)


def gust_strength(gusts, u):
    """
    Strength of a spot from one uniform draw u in [0, 1) (the same number of draws whatever the law, so that the
    other parts of the scenario do not change). Two laws:
      "strength": [low, high]          uniform between low and high (a spot can be weak);
      "strength_magnitude": [m1, m2]   a lull with probability "lull_probability" (0.5 by default), else a gust,
                                       of relative strength between m1 and m2 at its centre: every spot is clear.
    """
    if "strength_magnitude" not in gusts:
        low, high = map(float, gusts["strength"])
        return low + (high - low) * u
    m1, m2 = map(float, gusts["strength_magnitude"])
    p = float(gusts.get("lull_probability", 0.5))
    if u < p:
        return -(m1 + (m2 - m1) * u / p)
    return m1 + (m2 - m1) * (u - p) / (1.0 - p)


def max_neighbour_step(angle_deg):
    """Largest difference between two neighbouring cells (8-neighbourhood) of a smooth angle field, in degrees."""
    a = np.asarray(angle_deg, dtype=float)
    return float(max(np.abs(a[:, 1:] - a[:, :-1]).max(), np.abs(a[1:, :] - a[:-1, :]).max(),
                     np.abs(a[1:, 1:] - a[:-1, :-1]).max(), np.abs(a[1:, :-1] - a[:-1, 1:]).max()))


def wind_vectors(from_compass_deg, speed):
    """Wind vectors (blowing TO) from the compass direction the wind comes FROM (0 = north, clockwise)."""
    rad = np.radians(from_compass_deg)
    return np.stack([-speed * np.sin(rad), -speed * np.cos(rad)], axis=-1)


def rotate(field, delta_deg):
    """Turn every wind vector so that the direction it comes from moves clockwise by delta_deg (scalar or per cell)."""
    rad = np.radians(delta_deg)
    if np.ndim(rad) and field.dtype == np.float32:  # per-cell turn of a float32 field: float32 trigonometry is enough
        rad = rad.astype(np.float32)
    c, s = np.cos(rad), np.sin(rad)
    wx, wy = field[..., 0], field[..., 1]
    return np.stack([c * wx + s * wy, -s * wx + c * wy], axis=-1)


class WindModel:
    """Wind of one episode: field(t) is the wind field at step t, of shape (size, size, 2), indexed [y, x]."""

    def __init__(self, settings, seed, static=False, size=GRID_SIZE):
        self.settings, self.seed, self.static, self.size = settings, int(seed), bool(static), size
        self._make_front(settings.get("front"))
        self._make_base(settings["base"])
        self._make_oscillation(settings.get("oscillation"))
        self._make_gusts(settings.get("gusts"))
        turning = self.front is not None and self.front["kind"] == "turning"
        self.dynamic = (self.oscillating or self.gusting or turning) and not self.static
        self._last = None

    # ------------------------------------------------------------------ mechanism 1: base field
    def _make_base(self, base):
        seed, size = self.seed, self.size
        u = float(rng_stream(seed, STREAM_BASE_DIRECTION).random())
        lo, hi = map(float, base["mean_from_deg"])
        self.mean_from_deg = triangular(u, lo, hi) if base.get("mean_from_shape") == "triangular" else lo + (hi - lo) * u
        spatial = float(uniform(rng_stream(seed, STREAM_BASE_SPATIAL), *base["spatial_range_deg"]))
        small = smooth_field(rng_stream(seed, STREAM_BASE_SMALL), base["small_control"], base["small_std_deg"], size)
        limit = float(base["max_neighbour_deg"])
        front_now = 0.0
        if self.front is not None:  # the mean wind blows towards the side that looks good at the start
            f = self.front
            self.mean_from_deg = -f["side"] * f["mean_abs_deg"]
            spatial *= f["spatial_scale"]
            limit -= f["max_step"]
            front_now = self.front_twist(0)
        rng = rng_stream(seed, STREAM_BASE_LARGE)
        best = None
        for _ in range(MAX_PATTERN_DRAWS):
            g = unit_range_pattern(rng, base["spatial_control"], size)
            step = max_neighbour_step(spatial * g + small)
            if best is None or step < best[0]:
                best = (step, g)
            if step <= limit:
                break
        step, g = best
        while step > limit and spatial > 1.0:  # last resort: a weaker spatial variability
            spatial *= 0.95
            step = max_neighbour_step(spatial * g + small)
        self.spatial_range_deg = spatial
        self.base_from_deg = self.mean_from_deg + spatial * g + small + front_now
        self._base64 = wind_vectors(self.base_from_deg, float(base["speed"]))
        self._base = self._freeze(self._base64)

    # ------------------------------------------------------------------ mechanism 3: front
    def _make_front(self, front):
        """
        With probability `probability`, a front. side = +1 (right, x > 64) with probability `right_probability`:
        the mean wind then comes from the west (from -mean_abs_deg) and carries the boat to the right at the start.
        twist(x, y, t) = -a * tanh((x - 64.5) / width_x) * W(y, t), with two kinds of fronts:
          still    W = tanh((y0 - y) / width_y): good at the bottom of the lure side, bad higher up (y0 in start_y);
          turning  W = cos(pi * min(t / T, 1)): the whole lure side is good at the start and turns bad within T steps
                   (T in turning_steps); the other side does the opposite. Drawn with probability turning_probability.
        a, y0 and T are hidden.
        """
        self.front = None
        if front is None:
            return
        rng = rng_stream(self.seed, STREAM_FRONT)
        present, right, u_mean, u_twist, u_kind, u_y0, u_steps = (float(v) for v in rng.random(7))
        if present >= float(front["probability"]):
            return

        def draw(key, u):
            lo, hi = map(float, front[key])
            return lo + (hi - lo) * u
        turning = u_kind < float(front.get("turning_probability", 0.0))
        c = np.arange(self.size) + 0.5
        self.front = {"side": 1 if right < float(front["right_probability"]) else -1,
                      "kind": "turning" if turning else "still",
                      "mean_abs_deg": draw("mean_from_deg", u_mean), "twist_deg": draw("twist_deg", u_twist),
                      "y0": None if turning else draw("start_y", u_y0),
                      "turn_steps": draw("turning_steps", u_steps) if turning else None,
                      "width_y": float(front["width_y"]), "spatial_scale": float(front.get("spatial_scale", 1.0)),
                      "_sx": np.tanh((c - 0.5 * self.size) / float(front["width_x"])), "_y": c}
        # the steepest step of the twist alone, at any time: the base pattern gets the rest of the smoothness budget
        self.front["max_step"] = max_neighbour_step(self._twist(1.0 if turning else None))

    def _twist(self, w=None):
        """The twist for a uniform W = w (turning front), or for the still front line (w = None)."""
        f = self.front
        wy = np.full(self.size, w) if w is not None else np.tanh((f["y0"] - f["_y"]) / f["width_y"])
        return -f["twist_deg"] * wy[:, None] * f["_sx"][None, :]

    def front_twist(self, t):
        """Turn (degrees, per cell) that the front gives to the wind at step t; 0 without a front."""
        if self.front is None:
            return 0.0
        if self.front["kind"] == "still":
            return self._twist()
        return self._twist(math.cos(math.pi * min(float(t) / self.front["turn_steps"], 1.0)))

    # ------------------------------------------------------------------ mechanism 2: oscillation
    def _make_oscillation(self, osc):
        self.oscillating = osc is not None and not self.static
        if osc is None:
            self.period = self.amplitude_deg = self.phase = None
            return
        rng = rng_stream(self.seed, STREAM_OSCILLATION)
        self.period = float(uniform(rng, *osc["period_steps"]))
        self.amplitude_deg = float(uniform(rng, *osc["amplitude_deg"]))
        self.phase = float(uniform(rng, 0.0, 2.0 * math.pi))
        self._noise_std = float(osc["noise_deg"])
        self._noise_rho = float(osc["noise_correlation"])
        self._noise_rng = rng_stream(self.seed, STREAM_OSCILLATION_NOISE)
        self._noise = [self._noise_std * float(standard_normal(self._noise_rng, 1)[0])]

    def _noise_at(self, t):
        rho, std = self._noise_rho, self._noise_std
        while len(self._noise) <= t:  # always generated in the same order: deterministic
            z = float(standard_normal(self._noise_rng, 1)[0])
            self._noise.append(rho * self._noise[-1] + math.sqrt(1.0 - rho * rho) * std * z)
        return self._noise[t]

    def delta(self, t):
        """Rotation (degrees) of the whole field at step t with respect to the base field."""
        if not self.oscillating:
            return 0.0
        angle = 2.0 * math.pi * t / self.period + self.phase
        return self.amplitude_deg * math.sin(angle) + self._noise_at(int(t))

    # ------------------------------------------------------------------ mechanism 4: gusts and lulls
    def _make_gusts(self, gusts):
        self.gusting = gusts is not None and not self.static
        self.gusts = []
        if gusts is None:
            return
        rng = rng_stream(self.seed, STREAM_GUSTS)
        horizon, lead = float(gusts.get("horizon_steps", 700)), float(max(gusts["lifetime_steps"]))
        downwind = math.radians(self.mean_from_deg + 180.0)  # compass direction the mean wind blows to
        t = -lead
        while True:  # Poisson process of births, from uniform draws only
            t += -math.log(1.0 - float(rng.random())) / float(gusts["rate_per_step"])
            if t > horizon:
                break
            life = float(uniform(rng, *gusts["lifetime_steps"]))
            heading = downwind + math.radians(float(gusts["direction_spread_deg"]) * (2.0 * float(rng.random()) - 1.0))
            speed = float(uniform(rng, *gusts["speed_cells_per_step"]))
            self.gusts.append({
                "birth": t, "life": life,
                "x_mid": float(uniform(rng, 0.0, self.size)), "y_mid": float(uniform(rng, 0.0, self.size)),
                "vx": speed * math.sin(heading), "vy": speed * math.cos(heading),
                "radius": float(uniform(rng, *gusts["radius_cells"])),
                "strength": gust_strength(gusts, float(rng.random())),
                "shift": float(uniform(rng, *gusts.get("shift_deg", (0.0, 0.0)))),
            })
        self._gc = np.arange(self.size) + 0.5

    def gust_effect(self, t):
        """
        (relative speed factor, direction shift in degrees) per cell at step t. The spots combine multiplicatively:
        factor = prod (1 + strength) ** bump, so that the speed stays positive and varies smoothly in relative terms.
        A spot is a Gaussian bump, the outer product of two 1-D Gaussians: all the spots alive at step t are summed
        with one matrix product.
        """
        if not hasattr(self, "_g"):
            keys = ("birth", "life", "x_mid", "y_mid", "vx", "vy", "radius", "strength", "shift")
            self._g = {k: np.array([g[k] for g in self.gusts], dtype=float) for k in keys}
        g = self._g
        age = t - g["birth"]
        alive = (age >= 0.0) & (age <= g["life"])
        if not alive.any():
            return np.ones((self.size, self.size), dtype=np.float32), np.zeros((self.size, self.size), dtype=np.float32)
        envelope = np.sin(np.pi * age[alive] / g["life"][alive]) ** 2
        dt = (t - (g["birth"] + 0.5 * g["life"]))[alive]
        k = 1.0 / (2.0 * g["radius"][alive] ** 2)
        c = self._gc[None, :]
        ex = np.exp(-k[:, None] * (c - (g["x_mid"][alive] + g["vx"][alive] * dt)[:, None]) ** 2).astype(np.float32)
        ey = (envelope[:, None] * np.exp(-k[:, None] * (c - (g["y_mid"][alive] + g["vy"][alive] * dt)[:, None]) ** 2)
              ).astype(np.float32)
        log_factor = (ey * np.log1p(g["strength"][alive]).astype(np.float32)[:, None]).T @ ex
        shift = (ey * g["shift"][alive].astype(np.float32)[:, None]).T @ ex
        return np.exp(log_factor), shift

    # ------------------------------------------------------------------ fields
    def _turn_float32(self, turn):
        """rotate(base field, turn) for a per-cell turn, in float32 on the two components stored apart (faster)."""
        if not hasattr(self, "_bx"):
            self._bx, self._by = (np.ascontiguousarray(self._base[..., i]) for i in (0, 1))
        rad = np.asarray(turn, dtype=np.float32) * np.float32(math.pi / 180.0)
        c, s = np.cos(rad), np.sin(rad)
        out = np.empty(self._bx.shape + (2,), dtype=np.float32)
        out[..., 0] = c * self._bx + s * self._by
        out[..., 1] = c * self._by - s * self._bx
        return out

    @staticmethod
    def _freeze(field):
        field = np.ascontiguousarray(field, dtype=np.float32)
        field.setflags(write=False)
        return field

    def field(self, t):
        """Wind field at step t. Never modified after it is returned: an agent can keep past fields."""
        if not self.dynamic:
            return self._base
        t = int(t)
        if self._last is None or self._last[0] != t:
            turn = self.delta(t)
            if self.front is not None and self.front["kind"] == "turning":  # the base field holds the front at step 0
                if "_twist0" not in self.front:
                    self.front["_twist0"] = self.front_twist(0)
                turn = turn + (self.front_twist(t) - self.front["_twist0"])
            factor = None
            if self.gusting:
                factor, shift = self.gust_effect(t)
                turn = turn + shift
            if np.ndim(turn):  # per-cell turn: float32, the field is returned in float32
                f = self._turn_float32(turn)
            else:
                f = rotate(self._base64, turn)
            if factor is not None:
                f = f * factor[..., None]
            self._last = (t, self._freeze(f))
        return self._last[1]


def compass_from_deg(wind):
    """Compass direction the wind comes FROM, in degrees in [0, 360): 0 = from the north, 90 = from the east."""
    w = np.asarray(wind, dtype=float)
    return np.mod(np.degrees(np.arctan2(-w[..., 0], -w[..., 1])), 360.0)


def max_neighbour_differences(field):
    """Largest angle (degrees) and relative speed difference between two neighbouring cells (8-neighbourhood)."""
    w = np.asarray(field, dtype=float)
    pairs = [(w[:, :-1], w[:, 1:]),        # (y, x) and (y, x + 1)
             (w[:-1, :], w[1:, :]),        # (y, x) and (y + 1, x)
             (w[:-1, :-1], w[1:, 1:]),     # (y, x) and (y + 1, x + 1)
             (w[:-1, 1:], w[1:, :-1])]     # (y, x + 1) and (y + 1, x)
    worst_angle, worst_speed = 0.0, 0.0
    for a, b in pairs:
        sa, sb = np.hypot(a[..., 0], a[..., 1]), np.hypot(b[..., 0], b[..., 1])
        cos = np.clip((a * b).sum(-1) / np.maximum(sa * sb, 1e-12), -1.0, 1.0)
        worst_angle = max(worst_angle, float(np.degrees(np.arccos(cos)).max()))
        worst_speed = max(worst_speed, float((np.abs(sa - sb) / np.maximum(np.maximum(sa, sb), 1e-12)).max()))
    return worst_angle, worst_speed
