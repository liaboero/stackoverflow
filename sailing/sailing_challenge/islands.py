"""
Islands. An island is drawn from a smooth shape, but the rule is the set of forbidden cells:
a cell is forbidden when its centre lies inside the shape.

Layouts:
  central  one central island, same footprint as in v2 (symmetric about x = 64.5, the line through the start and the goal)
  three    a shorter central island, closer to the start, then one island on each side, of different sizes, in front
           of the inner side of the corridor exits: after choosing a side of the central island, choose again
  four     the same, plus a small island just south of the goal: the boat cannot always arrive from below
"""
import numpy as np

from .constants import GRID_SIZE

# Each island: a convex polygon (corners, in cells) whose corners are rounded with the given radii.
CENTRAL_ISLAND = {  # v2 footprint: tip at (64.5, 17), shoulders at y = 43, body over x in [38, 91) up to y = 86
    "corners": [(64.5, 17.0), (91.0, 43.0), (91.0, 86.0), (38.0, 86.0), (38.0, 43.0)],
    "radii": [4.0, 10.0, 8.0, 8.0, 10.0],
}
SHORT_CENTRAL_ISLAND = {  # shorter and closer to the start: tip at (64.5, 11), body over x in [40, 89) up to y = 58
    "corners": [(64.5, 11.0), (89.0, 28.0), (89.0, 58.0), (40.0, 58.0), (40.0, 28.0)],
    "radii": [3.0, 8.0, 7.0, 7.0, 8.0],
}
def scaled(island, factor):
    """The same shape, `factor` times larger around the mean of its corners."""
    corners = np.asarray(island["corners"], dtype=float)
    centre = corners.mean(axis=0)
    return {"corners": [tuple(c) for c in centre + factor * (corners - centre)],
            "radii": [factor * r for r in island["radii"]]}


LEFT_ISLAND = scaled({  # small, in front of the inner side of the left corridor exit: pass it on the inside or outside
    "corners": [(30.0, 72.0), (43.0, 74.0), (42.0, 88.0), (29.0, 86.0)],
    "radii": [4.0, 4.0, 4.0, 4.0],
}, 1.2)
RIGHT_ISLAND = scaled({  # larger and elongated, in front of the inner side of the right corridor exit
    "corners": [(82.0, 76.0), (104.0, 78.0), (105.0, 92.0), (83.0, 92.0)],
    "radii": [5.0, 5.0, 5.0, 5.0],
}, 1.2)
GOAL_ISLAND = scaled({  # small diamond just south of the goal: arrive from the left or from the right
    "corners": [(64.5, 105.0), (72.0, 112.0), (64.5, 119.0), (57.0, 112.0)],
    "radii": [3.0, 3.0, 3.0, 3.0],
}, 1.1)
LAYOUTS = {
    "central": {"islands": [CENTRAL_ISLAND], "symmetric": True},
    "three": {"islands": [SHORT_CENTRAL_ISLAND, LEFT_ISLAND, RIGHT_ISLAND], "symmetric": False},
    "four": {"islands": [SHORT_CENTRAL_ISLAND, LEFT_ISLAND, RIGHT_ISLAND, GOAL_ISLAND], "symmetric": False},
}


def rounded_polygon(corners, radii, points_per_arc=24):
    """Outline of a convex polygon whose corners are replaced by circle arcs of the given radii."""
    corners = np.asarray(corners, dtype=float)
    n = len(corners)
    outline = []
    for k in range(n):
        p, prev, nxt = corners[k], corners[k - 1], corners[(k + 1) % n]
        u1 = (prev - p) / np.linalg.norm(prev - p)
        u2 = (nxt - p) / np.linalg.norm(nxt - p)
        half = 0.5 * np.arccos(np.clip(np.dot(u1, u2), -1.0, 1.0))
        r = radii[k]
        bisector = (u1 + u2) / np.linalg.norm(u1 + u2)
        centre = p + bisector * r / np.sin(half)
        t1 = p + u1 * r / np.tan(half)  # tangent points on both edges
        t2 = p + u2 * r / np.tan(half)
        a1 = np.arctan2(*(t1 - centre)[::-1])
        a2 = np.arctan2(*(t2 - centre)[::-1])
        sweep = np.mod(a2 - a1 + np.pi, 2 * np.pi) - np.pi  # short arc from t1 to t2
        angles = a1 + sweep * np.linspace(0.0, 1.0, points_per_arc)
        outline.append(centre + r * np.stack([np.cos(angles), np.sin(angles)], axis=1))
    return np.concatenate(outline)


def points_in_polygon(px, py, polygon):
    """Even-odd rule, vectorised over the points."""
    inside = np.zeros(np.shape(px), dtype=bool)
    x, y = polygon[:, 0], polygon[:, 1]
    for k in range(len(polygon)):
        x1, y1, x2, y2 = x[k - 1], y[k - 1], x[k], y[k]
        if y1 == y2:
            continue
        crosses = (y1 > py) != (y2 > py)
        x_cross = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
        inside ^= crosses & (px < x_cross)
    return inside


def island_outlines(layout="central"):
    """Smooth outlines of the islands of a layout (for drawing). The forbidden cells are given by island_mask."""
    return [rounded_polygon(isl["corners"], isl["radii"]) for isl in LAYOUTS[layout]["islands"]]


def island_mask(layout="central", size=GRID_SIZE):
    """Forbidden cells, as a read-only boolean array indexed [y, x]."""
    centres = np.arange(size) + 0.5
    px, py = np.meshgrid(centres, centres)
    mask = np.zeros((size, size), dtype=bool)
    for outline in island_outlines(layout):
        mask |= points_in_polygon(px, py, outline)
    if LAYOUTS[layout]["symmetric"]:
        # exact left-right symmetry about x = 64.5: column j mirrors column 128 - j (column 0 has no mirror)
        mask[:, 1:size // 2] = mask[:, size - 1:size // 2:-1]
    mask.setflags(write=False)
    return mask
