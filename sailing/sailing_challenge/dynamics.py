"""
The exact one-step model, shared by the environment and by is_safe (a single implementation).

One step from position p and velocity v with an action, in the wind w of the boat's cell:
  1. v' = next_velocity(v, w, action); p' = p + v', kept inside the map: at a border the boat slides along it,
     and the component along the border (of v' and of the move) is multiplied by physics.wall_friction;
  2. if the segment [p, p'] enters the goal circle before any forbidden cell, the race is over: event "arrive",
     with `fraction` in (0, 1] the part of the step travelled before entering the circle;
  3. otherwise, if it enters a forbidden cell, the boat crashes: it stays at p with zero velocity (event "crash");
  4. otherwise it moves to p' with velocity v' (event "move").
"""
from typing import NamedTuple, Optional

import numpy as np

from .constants import GOAL, GOAL_RADIUS, GRID_SIZE
from .geometry import apply_borders, segment_circle_entry, segment_first_blocked
from .physics import Physics, next_velocity


class Transition(NamedTuple):
    event: str                 # "move", "crash" or "arrive"
    position: np.ndarray       # float64[2]
    velocity: np.ndarray       # float64[2]
    fraction: Optional[float]  # part of the step travelled before arriving (or before the crash), else None


def local_wind(wind_field, position):
    """Wind of the cell of `position` (float64[2]); wind_field is indexed [y, x]."""
    return np.asarray(wind_field[int(position[1]), int(position[0])], dtype=float)


def transition(position, velocity, wind, action, island_mask, physics=Physics(),
               goal=GOAL, goal_radius=GOAL_RADIUS, size=GRID_SIZE):
    p0 = np.asarray(position, dtype=float)
    v = next_velocity(velocity, wind, action, physics)
    px, py, vx, vy = apply_borders(p0[0], p0[1], v[0], v[1], size, physics.wall_friction)
    p1 = np.array([px, py])
    s_goal = segment_circle_entry(p0, p1, goal, goal_radius)
    s_hit = segment_first_blocked(p0, p1, island_mask)
    if s_goal is not None and (s_hit is None or s_goal <= s_hit):
        return Transition("arrive", p0 + s_goal * (p1 - p0), np.array([vx, vy]), s_goal)
    if s_hit is not None:
        return Transition("crash", p0.copy(), np.zeros(2), s_hit)
    return Transition("move", p1, np.array([vx, vy]), None)
