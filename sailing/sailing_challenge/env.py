"""
The sailing environment. Gymnasium interface when gymnasium is installed; the core only needs NumPy
(so that the same environment can also run in a web page).

    env = SailingEnv(regime="train")
    obs, info = env.reset(seed=3)            # info["env_info"]: fixed characteristics, given once per episode
    obs, reward, terminated, truncated, info = env.step(action)

Observation (a dict, at every step):
  position          float64[2]           exact continuous position (x east, y north)
  velocity          float64[2]           exact velocity
  local_wind        float32[2]           wind of the boat's cell (the one the next step uses)
  local_wind_speed  float32 (0-d array)  its speed
  local_wind_from   float32 (0-d array)  compass direction it comes from, degrees (0 = north, 90 = east)
  wind_field        float32[128,128,2]   wind of every cell, indexed [y, x] (a new array at every step)
  step              int                  number of steps played so far
The island mask is not repeated in the observation: it is in env_info.
"""
import secrets

import numpy as np

from .constants import (ACTION_NAMES, ACTION_VECTORS, DISCOUNT, GOAL, GOAL_RADIUS, GOAL_REWARD, GRID_SIZE,
                        MAX_STEPS, MAX_WIND_SPEED, N_ACTIONS, START, STAY, VERSION)
from .dynamics import local_wind, transition
from .islands import island_mask
from .physics import Physics
from .regimes import load_regime
from .wind import WindModel, compass_from_deg

try:
    import gymnasium
    from gymnasium import spaces
    _Base = gymnasium.Env
except ImportError:  # e.g. in the browser: same environment, without the gymnasium base class
    gymnasium, spaces, _Base = None, None, object


def parse_action(action):
    """(action, valid): an invalid action (None, not an integer, outside 0..8) is replaced by "stay"."""
    if isinstance(action, (bool, np.bool_)):
        return STAY, False
    if isinstance(action, (int, np.integer)) or (isinstance(action, np.ndarray) and action.ndim == 0
                                                 and np.issubdtype(action.dtype, np.integer)):
        a = int(action)
        if 0 <= a < N_ACTIONS:
            return a, True
    return STAY, False


class SailingEnv(_Base):
    metadata = {"render_modes": []}

    def __init__(self, regime="train", static_wind=False, physics=None, max_steps=MAX_STEPS):
        # regime="random": a new random regime (wind settings and island layout) is drawn at every reset(seed)
        self.random_regime = regime == "random"
        self.regime = load_regime("train" if self.random_regime else regime)
        self.static_wind = bool(static_wind)
        self.physics = physics if physics is not None else Physics()
        self.max_steps = int(max_steps)
        self.render_mode = None
        self.island_mask = island_mask(self.regime.get("islands", "central"))
        if spaces is not None:
            w_max = MAX_WIND_SPEED
            self.action_space = spaces.Discrete(N_ACTIONS)
            self.observation_space = spaces.Dict({
                "position": spaces.Box(0.0, float(GRID_SIZE), (2,), np.float64),
                "velocity": spaces.Box(-self.physics.max_speed, self.physics.max_speed, (2,), np.float64),
                "local_wind": spaces.Box(-w_max, w_max, (2,), np.float32),
                "local_wind_speed": spaces.Box(0.0, w_max, (), np.float32),
                "local_wind_from": spaces.Box(0.0, 360.0, (), np.float32),
                "wind_field": spaces.Box(-w_max, w_max, (GRID_SIZE, GRID_SIZE, 2), np.float32),
                "step": spaces.Discrete(self.max_steps + 1),
            })
        self.episode_seed = None
        self._done = True

    # ------------------------------------------------------------------ gymnasium API
    def reset(self, seed=None, options=None):
        if gymnasium is not None:
            super().reset(seed=seed)
        self.episode_seed = int(seed) if seed is not None else secrets.randbits(32)
        if self.random_regime:
            from .regimes.random_regime import random_regime
            self.regime = random_regime(self.episode_seed)
            self.island_mask = island_mask(self.regime.get("islands", "central"))
        self.wind = WindModel(self.regime["wind"], self.episode_seed, static=self.static_wind)
        self.t = 0
        self.position = np.array(START, dtype=float)
        self.velocity = np.zeros(2)
        self.n_crashes = 0
        self.n_invalid_actions = 0
        self.arrival_time = None
        self._field = self.wind.field(0)
        self._done = False
        return self._observation(), {"env_info": self.env_info}

    @property
    def env_info(self):
        """Fixed characteristics of the environment, given to the agent once per episode (a new dict at each call)."""
        return {
            "version": VERSION,
            "grid_size": GRID_SIZE,
            "start": START,
            "goal": GOAL,
            "goal_radius": GOAL_RADIUS,
            "island_mask": self.island_mask.copy(),  # forbidden cells, indexed [y, x]
            "action_names": ACTION_NAMES,
            "action_vectors": ACTION_VECTORS.copy(),
            "physics": self.physics.as_dict(),
            "max_steps": self.max_steps,
            "discount": DISCOUNT,
        }

    def step(self, action):
        if self._done:
            raise RuntimeError("the episode is over: call reset()")
        a, valid = parse_action(action)
        if not valid:
            self.n_invalid_actions += 1
        tr = transition(self.position, self.velocity, local_wind(self._field, self.position), a,
                        self.island_mask, self.physics)
        self.t += 1
        self.position, self.velocity = tr.position, tr.velocity
        reward, terminated = 0.0, False
        if tr.event == "arrive":
            reward, terminated = GOAL_REWARD, True
            self.arrival_time = (self.t - 1) + tr.fraction
        elif tr.event == "crash":
            self.n_crashes += 1
        truncated = not terminated and self.t >= self.max_steps
        self._done = terminated or truncated
        self._field = self.wind.field(self.t)
        info = {"event": tr.event, "crash": tr.event == "crash", "invalid_action": not valid,
                "arrival_time": self.arrival_time, "n_crashes": self.n_crashes,
                "n_invalid_actions": self.n_invalid_actions}
        return self._observation(), reward, terminated, truncated, info

    # ------------------------------------------------------------------ helpers
    def _observation(self):
        w = self._field[int(self.position[1]), int(self.position[0])]
        return {
            "position": self.position.copy(),
            "velocity": self.velocity.copy(),
            "local_wind": w.copy(),
            "local_wind_speed": np.array(np.hypot(w[0], w[1]), dtype=np.float32),
            "local_wind_from": np.array(compass_from_deg(w), dtype=np.float32),
            "wind_field": self._field.copy(),  # gymnasium: a new array at every call
            "step": self.t,
        }
