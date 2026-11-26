"""
Distance to the goal around the islands.

    distance = distance_to_goal(env_info)   # once: the map never changes
    distance[int(y), int(x)]                # at the boat's position (x, y), in cells

For every cell: the length of the shortest path from its centre to the goal, going from free cell to free cell in the
8 directions (1 for a straight move, sqrt(2) for a diagonal one), minus the goal radius. No rule decides on which side
of an island to go: from every cell, the shortest of all the paths is taken.
It knows nothing of the wind or of the boat: it tells where the boat can go through, not which route is the fastest.
A greedy agent that tries the 9 actions with next_step and keeps the one that lands closest to the goal by this
distance never gets stuck against an island (is_safe alone avoids the crashes, not the dead ends).
"""
import heapq
import math

import numpy as np

_MOVES = [(dx, dy, math.hypot(dx, dy)) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dx or dy]


def distance_to_goal(env_info):
    """
    Array of shape (grid_size, grid_size), indexed [y, x], in cells: 0 near the goal, inf on the islands.
    Takes a few hundredths of a second (Dijkstra's algorithm on the grid): compute it once, in reset() or __init__().
    """
    mask = np.asarray(env_info["island_mask"], dtype=bool)
    n = mask.shape[0]
    free = (~mask).ravel().tolist()
    gx, gy = (int(v) for v in env_info["goal"])
    dist = [math.inf] * (n * n)
    dist[gy * n + gx] = 0.0
    heap = [(0.0, gy * n + gx)]
    while heap:
        d, k = heapq.heappop(heap)
        if d > dist[k]:
            continue
        y, x = divmod(k, n)
        for dx, dy, cost in _MOVES:
            u, v = x + dx, y + dy
            if 0 <= u < n and 0 <= v < n:
                j = v * n + u
                if free[j] and d + cost < dist[j]:
                    dist[j] = d + cost
                    heapq.heappush(heap, (d + cost, j))
    return np.maximum(np.array(dist).reshape(n, n) - float(env_info["goal_radius"]), 0.0)
