# Sailing challenge: agents for wind maps never seen before

Test set: 75 episodes never used in any training. It holds the 5 game scenarios, 30 new seeds of the training
regime and 40 random regimes (wind settings and island layout drawn by `regimes/random_regime.py`). The table gives
the mean score, 100 · 0.995^steps.

| Agent | Scenarios | Training regime | Random regimes | All |
|---|---|---|---|---|
| ForecastPlanner, long training (64 episodes, `agent/trained_params.json`) | 72.41 | 73.47 | 71.57 | **72.39** |
| Same, deeper search at test time (horizon 7, beam 32) | 72.27 | 73.56 | 71.47 | 72.36 |
| ForecastPlanner, first training (32 episodes) | 72.65 | 73.48 | 70.96 | 72.08 |
| ForecastPlanner, hand-set parameters | 72.79 | 73.03 | 70.43 | 71.63 |
| Same planner, wind assumed constant | 70.63 | 70.49 | 69.99 | 70.24 |
| Greedy | 61.27 | 62.03 | 62.10 | 62.02 |
| AlphaStar-style network, imitation of the Oracle (`alphastar/`) | 66.69 | 66.32 | 63.38 | 64.77 |
| AlphaStar-style network, imitation then RL | 64.89 | 66.44 | 63.93 | 65.00 |
| "Oracle" reference (knows the future wind, limited beam search) | 73.20 | 72.69 | 71.78 | 72.24 |

The "Oracle" reference is not an upper bound: its search is limited (beam of 300). A wider search found 49.8 steps
on scenario 5, where this reference takes 59.8.

Reproduce: `python agent/train.py`, `python agent/evaluate.py`, `python alphastar/gen_demos.py --out demos.npz`,
`python alphastar/train_bc.py --demos demos.npz --out alphastar/bc.pt`,
`python alphastar/train_rl.py --init alphastar/bc.pt --out alphastar/rl.pt`,
`python agent/evaluate.py --agents Net:alphastar/rl.pt`.

## A* (hybrid A* over boat states, `agent/astar.py`)

A* was tested as an oracle, with the true future wind. It used the same Dijkstra time-to-goal map as heuristic and
merged states on a grid of 1 cell and 1 cell/step.

| Scenario | A*, weight 1 | A*, weight 1.2 | Tuned beam search |
|---|---|---|---|
| 1 | 54.4 | 53.6 | 49.7 |
| 2 | 80.8 | 85.0 | 77.8 |
| 3 | 63.6 | 64.5 | 58.9 |
| 4 | 61.6 | 61.6 | 61.4 |
| 5 | 57.4 | 61.2 | 49.8 |

A* needs 10,000 to 50,000 expansions per race and finds slower routes. The heuristic is not admissible (the wind
changes) and the state merging is coarse. With finer merging, A* stays stuck on head-wind plateaus.

## Weight of the turning-front forecast (`agent/sweep_front.py`)

Two trainings ended with front_weight on its upper bound (1.5). A third training allowed values up to 4. It found
nothing better than 74.61 in 4 generations, so it was stopped. A targeted test then varied only front_weight, on 28
episodes with a turning front and 28 without one, all new:

| front_weight | Turning front | Others |
|---|---|---|
| 1.0 | 75.63 | 69.84 |
| 1.5 (trained) | 75.19 | 69.69 |
| 2.0 | 73.45 | 69.71 |
| 2.5 | 73.40 | 69.81 |
| 3.0 | 71.08 | 69.65 |

Anticipating the front more strongly hurts (-1.7 ± 0.8 at 2.0, -4.1 ± 1.3 at 3.0). The bound was not what limited
the agent.
