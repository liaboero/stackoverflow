"""
Boat physics. Public, identical everywhere, and optional to use: an agent can also learn it by trial and error.

For an action of unit heading u and the wind w of the boat's cell:
  alpha  = angle between u and the direction the wind comes FROM, in [0, 180] degrees
  eta    = efficiency(alpha): 0.05 below 45 degrees (no-go zone), 0.5 -> 1 between 45 and 90,
           1 between 90 and 135, 1 -> 0.5 between 135 and 180
  v*     = eta * |w| * boat_performance * u          (target velocity; "stay" has v* = 0)
  v_new  = v* + inertia * (v - v*), then |v_new| <= max_speed
  p_new  = p + v_new                                   (no rounding)
All functions accept a single wind vector of shape (2,) or a whole field of shape (..., 2).
"""
from dataclasses import asdict, dataclass

import numpy as np

from .constants import ACTION_HEADINGS_DEG, ACTION_VECTORS, STAY


@dataclass(frozen=True)
class Physics:
    """Physical constants of the boat (public). The evaluation always uses the default values."""
    boat_performance: float = 0.4  # calibrated with the benchmark (20/11): a good race lasts about 50 to 100 steps
    inertia: float = 0.3
    max_speed: float = 8.0
    wall_friction: float = 0.5  # at a map border, the component along the border is multiplied by this factor

    def as_dict(self):
        return asdict(self)


def efficiency(alpha_deg):
    """Sailing efficiency for an angle alpha (degrees, in [0, 180]) between the heading and the wind-from direction."""
    a = np.asarray(alpha_deg, dtype=float)
    return np.where(a < 45.0, 0.05,
           np.where(a < 90.0, 0.5 + 0.5 * (a - 45.0) / 45.0,
           np.where(a < 135.0, 1.0,
                    np.maximum(0.5, 1.0 - 0.5 * (a - 135.0) / 45.0))))


def wind_from_deg(wind):
    """Direction the wind comes FROM, in degrees, mathematical convention (0 = from the east, 90 = from the north)."""
    w = np.asarray(wind, dtype=float)
    return np.degrees(np.arctan2(-w[..., 1], -w[..., 0]))


def angle_to_wind(wind, action):
    """Angle alpha in [0, 180] degrees between the heading of `action` (0..7) and the direction the wind comes from."""
    diff = np.mod(ACTION_HEADINGS_DEG[action] - wind_from_deg(wind) + 180.0, 360.0) - 180.0
    return np.abs(diff)


def target_velocity(wind, action, physics=Physics()):
    """Velocity v* the boat tends to when it keeps playing `action` in wind `wind`. Shape of `wind`: (..., 2)."""
    w = np.asarray(wind, dtype=float)
    if action == STAY:
        return np.zeros_like(w)
    speed = np.hypot(w[..., 0], w[..., 1])
    magnitude = efficiency(angle_to_wind(w, action)) * speed * physics.boat_performance
    return magnitude[..., None] * ACTION_VECTORS[action]


def target_velocities(wind, physics=Physics()):
    """Target velocities of the 9 actions at once: shape (9, ..., 2). Same values as target_velocity, action by action."""
    w = np.asarray(wind, dtype=float)
    speed = np.hypot(w[..., 0], w[..., 1])
    wind_from = wind_from_deg(w)
    out = np.zeros((len(ACTION_VECTORS),) + w.shape)
    for a in range(STAY):
        alpha = np.abs(np.mod(ACTION_HEADINGS_DEG[a] - wind_from + 180.0, 360.0) - 180.0)
        magnitude = efficiency(alpha) * speed * physics.boat_performance
        out[a] = magnitude[..., None] * ACTION_VECTORS[a]
    return out


def next_velocity(velocity, wind, action, physics=Physics()):
    """Velocity after one step from `velocity` with `action` in the local wind `wind` (inertia, then speed cap)."""
    v_target = target_velocity(wind, action, physics)
    v = v_target + physics.inertia * (np.asarray(velocity, dtype=float) - v_target)
    speed = np.hypot(v[..., 0], v[..., 1])
    scale = np.where(speed > physics.max_speed, physics.max_speed / np.maximum(speed, 1e-300), 1.0)
    return v * scale[..., None]
