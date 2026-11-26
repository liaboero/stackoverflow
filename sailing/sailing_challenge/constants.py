"""Fixed constants of the challenge. They are identical locally and on Codabench."""
import math

import numpy as np

VERSION = "3.0.0.dev0"

# The map is GRID_SIZE x GRID_SIZE cells of size 1. Positions are continuous, x points east and y points north.
# The cell of a position (x, y) is (floor(x), floor(y)); arrays over the map are indexed [y, x].
GRID_SIZE = 128
START = (64.5, 0.5)    # centre of cell (64, 0)
GOAL = (64.5, 127.5)   # centre of cell (64, 127)
GOAL_RADIUS = 2.0      # the race ends when the boat enters this circle

MAX_STEPS = 500
DISCOUNT = 0.995
GOAL_REWARD = 100.0
MAX_WIND_SPEED = 50.0  # bound of the observation space only (the base wind blows at 10)

# Actions: 8 headings and "stay". Headings are unit vectors, so diagonals are as fast as the other directions.
ACTION_NAMES = ("N", "NE", "E", "SE", "S", "SW", "W", "NW", "STAY")
N_ACTIONS = 9
STAY = 8
# Heading of actions 0..7 in degrees, mathematical convention (0 = east, 90 = north, counter-clockwise).
ACTION_HEADINGS_DEG = (90.0, 45.0, 0.0, -45.0, -90.0, -135.0, 180.0, 135.0)
_D = math.sqrt(0.5)
ACTION_VECTORS = np.array(
    [(0.0, 1.0), (_D, _D), (1.0, 0.0), (_D, -_D), (0.0, -1.0), (-_D, -_D), (-1.0, 0.0), (-_D, _D), (0.0, 0.0)])
ACTION_VECTORS.setflags(write=False)
