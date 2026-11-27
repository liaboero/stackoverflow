"""
A* over the states of the boat (hybrid A*): exact boat model, time-dependent wind.

A state is (position, velocity) at step t; its cost g is the number of steps played (t, and the fraction of the last
step for an arrival). The heuristic h is the time-to-goal map of the planner (Dijkstra on the grid from the goal),
read at the boat's position: f = g + weight * h. States are merged on a grid of positions and velocities (the
"hybrid" part: the search keeps exact continuous states, but visits each cell of that grid once).

With weight = 1 and an admissible map (one that never overestimates), A* would return the fastest route on the wind
it is given. Our map is an estimate (it ignores inertia, and the wind changes), so this is an efficient best-first
search rather than a proof of optimality; weight < 1 makes it more careful (more expansions), weight > 1 greedier.

    route = astar(position, velocity, t0, wind_at, mask, heuristic)   # wind_at(t) -> field [y, x, 2]
    route.actions, route.arrival_time, route.expansions
"""
import heapq
import itertools
from typing import List, NamedTuple, Optional

import numpy as np

from sailing_challenge.dynamics import transition


class Route(NamedTuple):
    actions: List[int]
    arrival_time: Optional[float]
    expansions: int


def astar(position, velocity, t0, wind_at, mask, heuristic, weight=1.0, max_expansions=30000,
          pos_res=0.5, vel_res=0.5, max_t=None):
    """heuristic(position, t) -> estimated steps to the goal."""
    counter = itertools.count()
    p0, v0 = np.asarray(position, float), np.asarray(velocity, float)
    start = (p0, v0, int(t0))
    heap = [(t0 + weight * heuristic(p0, t0), next(counter), start, None)]
    parents = {}     # node id -> (parent id, action)
    closed = set()
    best = None      # (arrival time, node id of the state before arriving, action)
    expansions = 0
    max_t = t0 + 400 if max_t is None else max_t
    while heap:
        f, nid, (p, v, t), _ = heapq.heappop(heap)
        if best is not None and f >= best[0]:
            break
        key = (round(p[0] / pos_res), round(p[1] / pos_res), round(v[0] / vel_res), round(v[1] / vel_res))
        if key in closed:
            continue
        closed.add(key)
        expansions += 1
        if expansions > max_expansions or t >= max_t:
            break
        field = wind_at(t)
        w = np.asarray(field[int(p[1]), int(p[0])], dtype=float)
        for a in range(9):
            tr = transition(p, v, w, a, mask)
            if tr.event == "arrive":
                at = t + tr.fraction
                if best is None or at < best[0]:
                    best = (at, nid, a)
            elif tr.event == "move":
                cid = next(counter)
                parents[cid] = (nid, a)
                heapq.heappush(heap, (t + 1 + weight * heuristic(tr.position, t + 1), cid,
                                      (tr.position, tr.velocity, t + 1), None))
    if best is None:
        # no arrival found within the budget: the route towards the most promising state that was expanded
        return Route([], None, expansions)
    actions = [best[2]]
    nid = best[1]
    while nid in parents:
        nid, a = parents[nid]
        actions.append(a)
    actions.reverse()
    return Route(actions, best[0], expansions)
