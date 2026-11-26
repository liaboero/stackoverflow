"""
Drawings with matplotlib, for the notebooks: the same look as the replay page (replay/index.html).

    from sailing_challenge.render import plot_env, plot_replay
    plot_env(env)                                  # the current state of an environment (after reset or step)
    plot_replay(replay, t=40)                      # a replay (build_replay / record_replay / load_replay) at step 40

pretty=True (default) paints the island (beach, forest, mountains, a lake) and a shaded sea; pretty=False shows the
plain map, where the forbidden cells are black. Either way the island covers exactly the forbidden cells: the pretty
island is outlined along the edges of these cells, beach included. Drawing only: nothing here changes the rules.

Needs matplotlib and scipy (not the environment itself, which only needs NumPy).
"""
import math

import numpy as np

from .constants import GOAL, GOAL_RADIUS, GRID_SIZE, START

TERRAIN_UPSAMPLE = 4
BOAT_COLORS = ("#d55e00", "#009e73", "#cc79a7", "#e69f00", "#7b3294", "#8c564b", "#b8860b", "#17a2b8")

# Colour ramps, identical to the replay page: (value, (r, g, b))
SEA_STOPS = ((0, (159, 224, 214)), (1.5, (108, 198, 207)), (5, (58, 156, 200)), (14, (36, 112, 168)),
             (30, (29, 90, 143)))
LAND_STOPS = ((0, (147, 194, 100)), (0.3, (106, 165, 72)), (0.55, (63, 122, 53)), (0.8, (141, 129, 115)),
              (1.0, (176, 167, 155)), (1.15, (244, 244, 242)))
SAND_WET, SAND_DRY = np.array((214, 189, 134.0)), np.array((236, 220, 174.0))
LAKE_DEEP, LAKE_SHORE = np.array((62, 134, 191.0)), np.array((121, 189, 227.0))
# Lake: off the centre of the deepest part of the island (offset in units of its depth), an elongated ellipse tilted
# by 25 degrees and a smaller lobe to the south-east, with an irregular shore. Same as the replay page.
LAKE_OFFSET = (0.18, 0.22)
LAKE_PARTS = ((0.0, 0.0, 1.35, 0.75, 25.0), (0.55, -0.75, 0.7, 0.55, -30.0))  # (dx/r, dy/r, a/r, b/r, angle in degrees)


# ---------------------------------------------------------------------------------------------- terrain (pretty)
def _hash2(ix, iy, seed):
    """Integer hash in [0, 1), the same as hash2 in the replay page (32-bit arithmetic)."""
    with np.errstate(over="ignore"):
        h = (np.asarray(ix).astype(np.uint32) * np.uint32(374761393)
             + np.asarray(iy).astype(np.uint32) * np.uint32(668265263)
             + np.uint32((seed * 982451653) % 2 ** 32))
        h = (h ^ (h >> np.uint32(13))) * np.uint32(1274126177)
        h = h ^ (h >> np.uint32(16))
    return h / 4294967296.0


def _value_noise(x, y, seed):
    x0, y0 = np.floor(x), np.floor(y)
    fx, fy = x - x0, y - y0
    u, v = fx * fx * (3 - 2 * fx), fy * fy * (3 - 2 * fy)
    x0, y0 = x0.astype(np.int64), y0.astype(np.int64)
    a, b = _hash2(x0, y0, seed), _hash2(x0 + 1, y0, seed)
    c, d = _hash2(x0, y0 + 1, seed), _hash2(x0 + 1, y0 + 1, seed)
    return a + (b - a) * u + (c - a) * v + (a - b - c + d) * u * v


def _fbm(x, y):
    return (0.5 * _value_noise(x / 8, y / 8, 1) + 0.3 * _value_noise(x / 4, y / 4, 2)
            + 0.2 * _value_noise(x / 2, y / 2, 3))


def _smoothstep(a, b, x):
    t = np.clip((x - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def _lake_shape(dx, dy, r):
    """1 on the shore of the lake (before the noise), < 1 inside: the smallest normalised distance to its parts."""
    e = np.full(np.shape(dx), np.inf)
    for px, py, a, b, angle in LAKE_PARTS:
        t = math.radians(angle)
        xx, yy = dx - px * r, dy - py * r
        u, v = xx * math.cos(t) + yy * math.sin(t), -xx * math.sin(t) + yy * math.cos(t)
        e = np.minimum(e, np.hypot(u / (a * r + 1e-9), v / (b * r + 1e-9)))
    return e


def _ramp(stops, x):
    xs = [s[0] for s in stops]
    return np.stack([np.interp(x, xs, [s[1][k] for s in stops]) for k in range(3)], axis=-1)


_terrain_cache = {}


def terrain_rgb(mask, upsample=TERRAIN_UPSAMPLE):
    """
    Pretty picture of the map, uint8 array (N, N, 3) with N = upsample * size, row 0 = south (use origin="lower").
    Every fine pixel of a forbidden cell is land (beach included), every other pixel is sea.
    """
    from scipy.ndimage import distance_transform_edt
    mask = np.asarray(mask, dtype=bool)
    key = (mask.tobytes(), mask.shape, upsample)
    if key in _terrain_cache:
        return _terrain_cache[key]
    F = upsample
    land = np.repeat(np.repeat(mask, F, axis=0), F, axis=1)
    N = land.shape[0]
    d_in = np.where(land, (distance_transform_edt(land) - 0.5) / F, 0.0)    # distance to the sea, in cells
    d_out = np.where(land, 0.0, (distance_transform_edt(~land) - 0.5) / F)  # distance to the island
    centres = (np.arange(N) + 0.5) / F
    x, y = np.meshgrid(centres, centres)
    noise = _fbm(x, y)
    dmax = float(d_in.max())
    deep = d_in >= 0.6 * dmax if dmax > 0 else np.zeros_like(land)
    cx, cy = ((float(x[deep].mean()) + LAKE_OFFSET[0] * dmax, float(y[deep].mean()) + LAKE_OFFSET[1] * dmax)
              if deep.any() else (0.0, 0.0))
    r = 0.2 * dmax
    e = _lake_shape(x - cx, y - cy, r) + 0.35 * (0.6 * _value_noise(x / 5, y / 5, 4)
                                                 + 0.4 * _value_noise(x / 2.5, y / 2.5, 5) - 0.5)
    lake = land & (r > 1) & (e < 1) & (d_in > 2.5)
    m = _smoothstep(0.0, 1.0, (e - 1) / 1.2)
    h = np.where(land, _smoothstep(0.0, 1.0, d_in / 14) * (0.3 + 0.7 * m) * (0.3 + 1.1 * noise * noise), 0.0)
    # sea
    sea = _ramp(SEA_STOPS, d_out) * (0.94 + 0.12 * noise)[..., None]
    foam = np.where(d_out < 0.3, (1 - d_out / 0.3) * 0.6, 0.0)[..., None]
    sea = sea + foam * (255 - sea)
    # land: sand on the edge, then a ramp on the elevation, with a hillshade (light from the north-west)
    gy, gx = np.gradient(h, 1.0 / F)
    shade = np.clip(1 - 3 * (gx * -math.sqrt(0.5) + gy * math.sqrt(0.5)), 0.55, 1.3)
    sand = (1 - _smoothstep(0.75, 1.25, d_in))[..., None]
    sand_c = SAND_WET + _smoothstep(0.0, 0.5, d_in)[..., None] * (SAND_DRY - SAND_WET)
    land_c = (sand * sand_c + (1 - sand) * _ramp(LAND_STOPS, h)) * np.where(sand > 0.5, 1.0, shade[..., None])
    lake_c = LAKE_DEEP + _smoothstep(0.6, 1.0, e)[..., None] * (LAKE_SHORE - LAKE_DEEP)
    rgb = np.where(land[..., None], np.where(lake[..., None], lake_c, land_c), sea)
    rgb = np.clip(np.rint(rgb), 0, 255).astype(np.uint8)
    rgb.setflags(write=False)
    _terrain_cache.clear()  # one map at a time is enough
    _terrain_cache[key] = rgb
    return rgb


def coast_segments(mask):
    """Edges between a forbidden cell and a free one (or the border), as a list of ((x0, y0), (x1, y1))."""
    mask = np.asarray(mask, dtype=bool)
    p = np.pad(mask, 1)
    segs = []
    for j, i in zip(*np.nonzero(mask)):
        jj, ii = j + 1, i + 1
        if not p[jj, ii - 1]:
            segs.append(((i, j), (i, j + 1)))
        if not p[jj, ii + 1]:
            segs.append(((i + 1, j), (i + 1, j + 1)))
        if not p[jj - 1, ii]:
            segs.append(((i, j), (i + 1, j)))
        if not p[jj + 1, ii]:
            segs.append(((i, j + 1), (i + 1, j + 1)))
    return segs


# ---------------------------------------------------------------------------------------------- drawing
def _downsample(field, resolution):
    from .replay import downsample_wind
    return downsample_wind(field, resolution)


def _new_axes(ax, size=6.0):
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(size, size))
    return ax


GOAL_RING = (0.90, 0.52, 0.06, 0.95)   # the edge of the goal circle (pretty mode)
GLOW_EXTENT = 5.0                      # the glow reaches this many goal radii


def goal_glow_rgba(n=128):
    """The light of the lighthouse: an RGBA image (n, n, 4) of the square GLOW_EXTENT goal radii around the goal.

    Concentric: pale yellow at the centre, orange at the edge of the circle, a soft halo beyond it."""
    c = (np.arange(n) + 0.5) / n * 2 - 1
    r = np.hypot(*np.meshgrid(c, c)) * GLOW_EXTENT            # in goal radii
    t = np.clip(r, 0, 1)[..., None]
    centre, edge = np.array([1.0, 0.97, 0.74]), np.array([1.0, 0.74, 0.26])
    rgba = np.zeros((n, n, 4))
    rgba[..., :3] = centre * (1 - t) + edge * t
    halo = 0.80 * np.exp(-(np.maximum(r - 1, 0) / 2.6) ** 1.2) * np.clip(1 - (r / GLOW_EXTENT) ** 3, 0, 1)
    rgba[..., 3] = np.where(r <= 1, 0.92 - 0.20 * r, halo)
    return rgba


def _glow_extent():
    e = GLOW_EXTENT * GOAL_RADIUS
    return (GOAL[0] - e, GOAL[0] + e, GOAL[1] - e, GOAL[1] + e)


def plot_map(mask, wind=None, ax=None, pretty=True, arrows=True, size=GRID_SIZE):
    """
    The map: sea, island, goal, start, and the wind if given.
    wind: a field (size, size, 2) indexed [y, x], or an already downsampled one (n, n, 2).
    """
    import matplotlib.colors as mcolors
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Circle, Rectangle
    ax = _new_axes(ax)
    mask = np.asarray(mask, dtype=bool)
    extent = (0, size, 0, size)
    small = None
    if wind is not None:
        wind = np.asarray(wind, dtype=float)
        small = _downsample(wind, 32) if wind.shape[0] == size and size % 32 == 0 else wind
    speed = None if small is None else np.hypot(small[..., 0], small[..., 1])
    if pretty:
        ax.imshow(terrain_rgb(mask), origin="lower", extent=extent, interpolation="bilinear", zorder=0)
        if speed is not None:  # whitecaps where the wind is stronger than 10, darker water where it is lighter
            rgba = np.zeros(speed.shape + (4,))
            strong = speed >= 10
            rgba[strong] = (1, 1, 1, 0)
            rgba[~strong] = (10 / 255, 30 / 255, 60 / 255, 0)
            rgba[..., 3] = np.where(strong, np.minimum(0.45, (speed - 10) / 10 * 0.45),
                                    np.minimum(0.4, (10 - speed) / 10 * 0.4))
            ax.imshow(rgba, origin="lower", extent=extent, interpolation="bilinear", zorder=1)
        ax.add_collection(LineCollection(coast_segments(mask), colors=[(70 / 255, 50 / 255, 25 / 255, 0.9)],
                                         linewidths=0.8, zorder=2))
    else:
        if speed is not None:
            cmap = mcolors.LinearSegmentedColormap.from_list("wind", ["#eef5fb", "#6f9fcb"])
            ax.imshow(speed, origin="lower", extent=extent, cmap=cmap, vmin=0, vmax=20, interpolation="bilinear",
                      zorder=0)
        else:
            ax.set_facecolor("#eef5fb")
        black = np.zeros(mask.shape + (4,))
        black[mask, 3] = 1.0
        ax.imshow(black, origin="lower", extent=extent, interpolation="nearest", zorder=1)
    if small is not None and arrows:
        n = small.shape[0]
        b = size / n
        stride = 2 if n > 16 else 1
        c = (np.arange(n) + 0.5) * b
        X, Y = np.meshgrid(c, c)
        U, V = small[..., 0].copy(), small[..., 1].copy()
        keep = ~mask[np.floor(Y).astype(int), np.floor(X).astype(int)]
        sel = (slice(stride // 2, None, stride),) * 2
        k = keep[sel]
        ax.quiver(X[sel][k], Y[sel][k], U[sel][k], V[sel][k], color=(1, 1, 1, 0.65) if pretty else "#193452",
                  angles="xy", scale_units="xy", scale=10 / (0.75 * b * stride), width=0.003, pivot="middle",
                  zorder=3)
    if pretty:  # the goal as the light of a lighthouse: an orange-yellow glow, brighter at the centre
        ax.imshow(goal_glow_rgba(), origin="lower", extent=_glow_extent(), interpolation="bilinear", zorder=4)
    ax.add_patch(Circle(GOAL, GOAL_RADIUS, facecolor="none", edgecolor=GOAL_RING if pretty else "#c99700",
                        linewidth=1.2 if pretty else 2, zorder=4))
    ax.add_patch(Rectangle((START[0] - 0.9, START[1] - 0.5), 1.8, 1.0, color="#7a5a3a" if pretty else "#3b4654",
                           zorder=4))
    ax.set_xlim(0, size)
    ax.set_ylim(0, size)
    ax.set_aspect("equal")
    ax.set_xticks(range(0, size + 1, 32))
    ax.set_yticks(range(0, size + 1, 32))
    return ax


def draw_boat(ax, position, heading, color, wind=(0.0, -1.0), pretty=True, length=4.5, alpha=1.0):
    """A boat at `position`, pointing along `heading`. Pretty: hull, white sail on the leeward side, wake."""
    from matplotlib.patches import Polygon
    x, y = map(float, position)
    hx, hy = map(float, heading)
    hn = math.hypot(hx, hy)
    moving = hn > 0.05
    hx, hy = (0.0, 1.0) if hn < 1e-9 else (hx / hn, hy / hn)
    L = length

    def P(a, c):  # a: along the hull, c: to starboard
        return (x + hx * a + hy * c, y + hy * a - hx * c)

    if not pretty:
        W = 0.5 * L * 0.6
        ax.add_patch(Polygon([P(0.6 * L * 0.6, 0), P(-0.4 * L * 0.6, W), P(-0.4 * L * 0.6, -W)], closed=True,
                             facecolor=color, edgecolor="white", linewidth=0.8, alpha=alpha, zorder=6))
        return
    W = 0.36 * L
    if moving:
        for side in (1, -1):
            (x0, y0), (x1, y1) = P(-0.42 * L, side * 0.12 * L), P(-1.0 * L, side * 0.28 * L)
            ax.plot([x0, x1], [y0, y1], color=(1, 1, 1, 0.55 * alpha), linewidth=1.2, zorder=5)
    t = np.linspace(0, 1, 12)[:, None]
    bow, c1, e1 = np.array(P(0.6 * L, 0)), np.array(P(0.25 * L, W * 0.62)), np.array(P(-0.4 * L, W / 2))
    c2, e2 = np.array(P(0.25 * L, -W * 0.62)), np.array(P(-0.4 * L, -W / 2))
    side1 = (1 - t) ** 2 * bow + 2 * (1 - t) * t * c1 + t ** 2 * e1
    side2 = (1 - t) ** 2 * e2 + 2 * (1 - t) * t * c2 + t ** 2 * bow
    hull = np.vstack([side1, P(-0.45 * L, 0.38 * W), P(-0.45 * L, -0.38 * W), side2])
    ax.add_patch(Polygon(hull, closed=True, facecolor=color, edgecolor=(0.12, 0.08, 0.04, 0.8), linewidth=0.7,
                         alpha=alpha, zorder=6))
    wn = math.hypot(*wind) or 1.0
    wx, wy = wind[0] / wn, wind[1] / wn
    side = -1 if (hx * wy - hy * wx) > 0 else 1
    along = hx * wx + hy * wy
    angle = (0.35 + 0.55 * (along + 1) / 2) * side
    mast = np.array(P(0.15 * L, 0))
    boom_len = 0.72 * L
    ca, sa = math.cos(angle), math.sin(angle)
    bx, by = -hx * ca + hy * sa, -hy * ca - hx * sa
    clew = mast + boom_len * np.array((bx, by))
    nx, ny = by, -bx
    if nx * wx + ny * wy < 0:
        nx, ny = -nx, -ny
    belly = mast + boom_len * 0.45 * np.array((bx, by)) + 0.3 * L * np.array((nx, ny))
    sail = (1 - t) ** 2 * mast + 2 * (1 - t) * t * belly + t ** 2 * clew
    ax.add_patch(Polygon(sail, closed=True, facecolor=(1, 1, 1, 0.97 * alpha), edgecolor=(0.16, 0.16, 0.16, 0.9),
                         linewidth=0.6, zorder=7))
    ax.plot(*mast, "o", color="#3b2a1a", markersize=1.5, alpha=alpha, zorder=8)


def draw_step_counter(ax, step):
    """The current step, large, in the top-left corner of the map (as on the replay page)."""
    ax.text(0.02, 0.98, f"Step {int(math.floor(float(step) + 1e-9))}", transform=ax.transAxes, ha="left", va="top",
            fontsize=15, fontweight="bold", color="white", zorder=10,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": (10 / 255, 20 / 255, 35 / 255, 0.6), "edgecolor": "none"})


def _wind_at(small, size, position):
    n = small.shape[0]
    b = size / n
    i = min(n - 1, max(0, int(position[0] // b)))
    j = min(n - 1, max(0, int(position[1] // b)))
    return tuple(small[j, i])


def plot_env(env, trajectory=None, ax=None, pretty=True, arrows=True, color=BOAT_COLORS[0]):
    """The current state of an environment: map, wind of the current step, boat, and a trajectory if given."""
    ax = plot_map(env.island_mask, env._field, ax=ax, pretty=pretty, arrows=arrows)
    small = _downsample(env._field, 32)
    if trajectory is not None and len(trajectory) > 1:
        tr = np.asarray(trajectory)
        ax.plot(tr[:, 0], tr[:, 1], color=color, linewidth=1.8, alpha=0.9, zorder=5)
    heading = env.velocity if np.hypot(*env.velocity) > 1e-9 else (0.0, 1.0)
    draw_boat(ax, env.position, heading, color, _wind_at(small, GRID_SIZE, env.position), pretty)
    ax.set_title(f"step {env.t}")
    draw_step_counter(ax, env.t)
    return ax


def plot_replay(replay, t=None, ax=None, pretty=True, arrows=True, legend=True):
    """A replay at step t (default: the end of the longest race): map, wind of that step, trails and boats."""
    from .replay import decode_mask, replay_wind_at
    ax = _new_axes(ax)
    size = replay["grid_size"]
    boats = replay["boats"]
    t_end = max(b["arrival_time"] if b["arrived"] else b["steps"] for b in boats)
    t = t_end if t is None else min(max(float(t), 0.0), t_end)
    small = replay_wind_at(replay, int(math.floor(t + 1e-9)))
    plot_map(decode_mask(replay["island_mask"], size), small, ax=ax, pretty=pretty, arrows=arrows, size=size)
    for k, b in enumerate(boats):
        color = BOAT_COLORS[k % len(BOAT_COLORS)]
        tr = np.asarray(b["trajectory"], dtype=float)
        times = np.arange(len(tr), dtype=float)
        if b["arrived"] and b["arrival_time"] is not None and len(tr) > 1:
            times[-1] = min(times[-1], max(times[-2], b["arrival_time"]))
        if t >= times[-1] or len(tr) == 1:
            pos = tr[-1]
            heading = tr[-1] - tr[-2] if len(tr) > 1 else np.array((0.0, 1.0))
            trail = tr
        else:
            k0 = int(np.searchsorted(times, t, side="right") - 1)
            f = (t - times[k0]) / (times[k0 + 1] - times[k0])
            pos = tr[k0] + f * (tr[k0 + 1] - tr[k0])
            heading = tr[k0 + 1] - tr[k0]
            trail = np.vstack([tr[:k0 + 1], pos])
        label = b["name"] + (f"  {b['score']:.1f}" if t >= times[-1] else "")
        ax.plot(trail[:, 0], trail[:, 1], color=color, linewidth=1.8, alpha=0.9, zorder=5, label=label)
        draw_boat(ax, pos, heading, color, _wind_at(small, size, pos), pretty)
    if legend:
        ax.legend(loc="lower right", fontsize=8, framealpha=0.85)  # the step counter is in the top-left corner
    ax.set_title(f"seed {replay['seed']} · step {t:.1f}")
    draw_step_counter(ax, t)
    return ax
