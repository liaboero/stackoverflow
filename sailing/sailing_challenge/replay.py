"""
Replays: export episodes to a compact JSON file for the replay page (replay/index.html).

The simplest way, from a notebook or a script (one file per seed, all the agents race on it):

    from sailing_challenge import export_replays
    export_replays({"greedy": GreedyAgent(), "mine": MyAgent()}, seeds=[1, 2, 3], out_dir="replays")

and from a terminal, with submission folders (agent.py with a class MyAgent):

    python -m sailing_challenge.make_replays my_submission/ other_submission/ --seeds 1-3 --out replays

Then open replay/index.html in a browser and choose the files. A single seed:

    record_replay({"greedy": GreedyAgent(), "mine": MyAgent()}, seed=3, path="replay_seed3.json")

or, from results of run_episode(..., record=True) on the same seed (one per boat):

    results = [run_episode(agent, seed=3, record=True) for agent in agents]
    save_replay(build_replay(results, names=["greedy", "mine"]), "replay_seed3.json")

The wind is exogenous (it does not depend on the boats), so it is rebuilt by replaying the environment
on the same seed: every boat of a replay sails in the same wind.

Replays are meant for the training regime: by default, build_replay refuses any other regime.

Format (version 1), one episode per file:
  format, version, package_version     "sailing-replay", 1, version of the package
  regime, seed, static_wind
  grid_size, start, goal, goal_radius, max_steps, physics
  island_mask    run lengths of the flattened [y, x] mask, starting with a run of free cells
  wind           {"resolution": n, "scale": s, "frames": [{"t": t, "u": [...], "v": [...]}, ...]}
                 the wind field averaged over blocks of (grid_size / n)^2 cells, as integers (value * s),
                 rows from south to north, west to east (index [j * n + i] for block row j, column i).
                 A frame holds the field from step t until the next frame (a frame identical to the
                 previous one is not repeated): the field observed at step t, used by the step t -> t + 1.
  boats          [{"name", "score", "arrived", "arrival_time", "steps", "crashes", "invalid_actions",
                   "trajectory": [[x, y], ...]}]   trajectory[k] is the position after k steps.
"""
import copy
import json
import os

import numpy as np

from .constants import GOAL, GOAL_RADIUS, GRID_SIZE, START, STAY, VERSION
from .env import SailingEnv
from .regimes import load_regime

FORMAT = "sailing-replay"
FORMAT_VERSION = 1
POSITION_DECIMALS = 3
REPLAY_PAGE_URL = "https://trahier.github.io/sailing-challenge/replay/"  # GitHub Pages of the public repository


# ---------------------------------------------------------------------------------------------- encoding helpers
def encode_mask(mask):
    """Run lengths of a boolean array (flattened row-major), starting with a run of False (possibly empty)."""
    flat = np.asarray(mask, dtype=bool).ravel()
    change = np.flatnonzero(flat[1:] != flat[:-1]) + 1
    bounds = np.concatenate([[0], change, [flat.size]])
    runs = np.diff(bounds).tolist()
    if flat.size and flat[0]:
        runs = [0] + runs
    return runs


def decode_mask(runs, size=GRID_SIZE):
    """Inverse of encode_mask: boolean array of shape (size, size)."""
    values = np.arange(len(runs)) % 2 == 1
    flat = np.repeat(values, runs)
    if flat.size != size * size:
        raise ValueError(f"the mask has {flat.size} cells, expected {size * size}")
    return flat.reshape(size, size)


def downsample_wind(field, resolution=32):
    """Mean of the wind vectors over square blocks: (size, size, 2) -> (resolution, resolution, 2), indexed [y, x]."""
    field = np.asarray(field, dtype=float)
    size = field.shape[0]
    if size % resolution:
        raise ValueError(f"the wind resolution ({resolution}) must divide the grid size ({size})")
    b = size // resolution
    return field.reshape(resolution, b, resolution, b, 2).mean(axis=(1, 3))


def _quantize(field, resolution, scale):
    small = downsample_wind(field, resolution)
    return np.rint(small * scale).astype(np.int32)


def wind_frames(fields, resolution=32, decimals=1):
    """
    Frames of the wind section from the fields at steps 0, 1, 2... (an iterable of (size, size, 2) arrays).
    Each frame is the field rounded to `decimals` and averaged over blocks; identical consecutive frames are dropped.
    """
    scale = 10 ** decimals
    frames, previous = [], None
    for t, field in enumerate(fields):
        q = _quantize(field, resolution, scale)
        if previous is not None and np.array_equal(q, previous):
            continue
        frames.append({"t": t, "u": q[..., 0].ravel().tolist(), "v": q[..., 1].ravel().tolist()})
        previous = q
    return {"resolution": resolution, "scale": scale, "frames": frames}


def episode_wind_fields(seed, n_steps, regime="train", static_wind=False, physics=None):
    """
    Wind fields at steps 0..n_steps of an episode, obtained by replaying the environment with "stay":
    the wind does not depend on the actions, so they are the fields every boat met on this seed.
    """
    env = SailingEnv(regime, static_wind=static_wind, physics=physics)
    obs, _ = env.reset(seed=seed)
    yield obs["wind_field"]
    for _ in range(min(int(n_steps), env.max_steps)):
        obs, _, terminated, truncated, _ = env.step(STAY)
        yield obs["wind_field"]
        if terminated:  # cannot happen with "stay" from the start, kept for safety
            break


# ---------------------------------------------------------------------------------------------- building a replay
def _check_regime(regime, allow_other_regimes):
    if regime == "train" or allow_other_regimes:
        return
    raise ValueError("replays are meant for the training regime only (regime='train'); "
                     "pass allow_other_regimes=True for a regime of your own")


def build_replay(results, names=None, regime="train", static_wind=False, physics=None, wind_resolution=32,
                 wind_decimals=1, allow_other_regimes=False):
    """
    Replay (a JSON-serialisable dict) from results of run_episode(..., record=True) on the same seed.

    results         list of results, or dict {name: result}
    names           names of the boats (default: "boat 1", "boat 2"...), ignored when results is a dict
    regime, static_wind, physics   the settings the episodes were run with (used to rebuild the wind)
    wind_resolution blocks per side of the wind field in the file (must divide 128)
    wind_decimals   rounding of the wind vectors
    """
    _check_regime(regime, allow_other_regimes)
    if isinstance(results, dict):
        names, results = list(results), list(results.values())
    results = list(results)
    if not results:
        raise ValueError("no result to replay")
    if names is None:
        names = [f"boat {k + 1}" for k in range(len(results))]
    if len(names) != len(results):
        raise ValueError("one name per result is needed")
    seeds = {r["seed"] for r in results}
    if len(seeds) != 1:
        raise ValueError(f"all the results must come from the same seed, got {sorted(seeds)}")
    seed = seeds.pop()
    if seed is None:
        raise ValueError("the seed of the episode is needed to rebuild the wind")
    boats = []
    for name, r in zip(names, results):
        if not r.get("trajectory"):
            raise ValueError(f"no trajectory for {name!r}: run the episode with record=True")
        traj = np.round(np.asarray(r["trajectory"], dtype=float), POSITION_DECIMALS)
        boats.append({
            "name": str(name),
            "score": round(float(r["score"]), 4),
            "arrived": bool(r["arrived"]),
            "arrival_time": None if r["arrival_time"] is None else round(float(r["arrival_time"]), 4),
            "steps": int(r["steps"]),
            "crashes": int(r["crashes"]),
            "invalid_actions": int(r["invalid_actions"]),
            "trajectory": traj.tolist(),
        })
    n_steps = max(b["steps"] for b in boats)
    env = SailingEnv(regime, static_wind=static_wind, physics=physics)
    wind = wind_frames(episode_wind_fields(seed, n_steps, regime, static_wind, physics),
                       wind_resolution, wind_decimals)
    regime_name = regime if isinstance(regime, str) else env.regime.get("name", "custom")
    return {
        "format": FORMAT,
        "version": FORMAT_VERSION,
        "package_version": VERSION,
        "regime": regime_name,
        "seed": int(seed),
        "static_wind": bool(static_wind),
        "grid_size": GRID_SIZE,
        "start": list(START),
        "goal": list(GOAL),
        "goal_radius": GOAL_RADIUS,
        "max_steps": env.max_steps,
        "physics": env.physics.as_dict(),
        "island_mask": encode_mask(env.island_mask),
        "wind": wind,
        "boats": boats,
    }


def save_replay(replay, path):
    """Write a replay as compact JSON."""
    with open(path, "w") as f:
        json.dump(replay, f, separators=(",", ":"))


def load_replay(path):
    with open(path) as f:
        replay = json.load(f)
    if replay.get("format") != FORMAT:
        raise ValueError(f"{path} is not a sailing replay")
    return replay


def record_replay(agents, seed, path=None, regime="train", static_wind=False, physics=None, **kwargs):
    """
    Run every agent of `agents` ({name: agent}) on the same seed and build the replay (written to `path` if given).
    Extra keyword arguments go to build_replay.
    """
    from .evaluation import run_episode
    results = {name: run_episode(agent, seed, regime=regime, static_wind=static_wind, physics=physics, record=True)
               for name, agent in agents.items()}
    replay = build_replay(results, regime=regime, static_wind=static_wind, physics=physics, **kwargs)
    if path is not None:
        save_replay(replay, path)
    return replay


def replay_page():
    """Where the replay page is: the file replay/index.html of the repository, else its web address."""
    local = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "replay", "index.html")
    return local if os.path.isfile(local) else REPLAY_PAGE_URL


def export_replays(agents, seeds=None, out_dir="replays", regime="train", static_wind=False, verbose=True, **kwargs):
    """
    Race the agents on each seed and write one replay file per seed, out_dir/replay_seed<seed>.json.
    Open the replay page (see replay_page()) and choose these files to watch the races.

    agents   {name: agent}, or a single agent (named "agent")
    seeds    the seeds to race on; default: the first 5 fixed seeds of the training regime
    Every episode plays a fresh copy of the agent (copy.deepcopy), as in the evaluation.
    Extra keyword arguments go to build_replay. Returns the paths of the files.
    """
    if not isinstance(agents, dict):
        agents = {"agent": agents}
    if seeds is None:
        seeds = load_regime(regime)["fixed_seeds"][:5]
    from .evaluation import run_episode
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for seed in seeds:
        results = {name: run_episode(copy.deepcopy(agent), int(seed), regime=regime, static_wind=static_wind,
                                     record=True)
                   for name, agent in agents.items()}
        replay = build_replay(results, regime=regime, static_wind=static_wind, **kwargs)
        path = os.path.join(out_dir, f"replay_seed{int(seed)}.json")
        save_replay(replay, path)
        paths.append(path)
        if verbose:
            scores = ", ".join(f"{b['name']} {b['score']:.1f}" for b in replay["boats"])
            print(f"seed {int(seed)}: {scores} -> {path}")
    if verbose:
        print(f"\nOpen the replay page in a browser and choose these files:\n  {replay_page()}")
    return paths


def replay_wind_at(replay, t):
    """Wind field of the replay at step t, as an array (resolution, resolution, 2) indexed [y, x] (for checks)."""
    wind = replay["wind"]
    frame = [f for f in wind["frames"] if f["t"] <= t][-1]
    n = wind["resolution"]
    u = np.asarray(frame["u"], dtype=float).reshape(n, n)
    v = np.asarray(frame["v"], dtype=float).reshape(n, n)
    return np.stack([u, v], axis=-1) / wind["scale"]
