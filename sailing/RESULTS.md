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

## Improving the forecast, step by step (`agent/diagnose.py`, `agent/sweep.py`)

**1. Diagnostic.** The planner was given the true future of one wind component, with its own forecast for the
others. Gain over the agent, on 48 development episodes:

| True future of... | Gain |
|---|---|
| global rotation (oscillation) | +1.16 ± 0.48 |
| gusts (ignored until now) | +1.12 ± 0.56 |
| turning front | +0.38 ± 0.20 |
| everything | +1.69 ± 0.57 |

**2. Gust forecast.** The drift of the speed pattern is measured by a global Lucas-Kanade fit between two steps and
extrapolated, with gusts fading over 30 steps. A first version filled the area entering from the map edge with a copy
of the edge; a long shift then spread one gust over a whole band of the map (scenario 2: 83 → 102 steps). It now fills
that area with "no gust".

**3. Oscillation forecast.** A Bayesian average over 7,488 sinusoids (period 25-140, amplitude 0-60, phase), weighted
by their fit to the measured rotation, replaces the single sinusoid fitted after 19 steps. Error 10 steps ahead, at
step 10: 5.7° instead of 13.2°. Most of the gain the diagnostic promised is out of reach, though: 30 steps ahead the
rotation cannot be known early in the race (error 25° at step 10, for every method). A prior narrowed to the
training regime was worse (-0.2).

**4. Hedging on the front** (a mix of the maps "the front turns" and "it does not"): no gain (-0.03 to -1.9). Not kept.

**5. Time-dependent map** (every cell uses the forecast at the time the boat would pass there): worse (-0.4 to
-0.9). Not kept.

**Test set (75 episodes), paired against the previous best (72.39):**

| Version | Scenarios | Training regime | Random regimes | All | Gain |
|---|---|---|---|---|---|
| A: new oscillation forecast | 72.51 | 73.33 | 72.22 | 72.68 | +0.30 ± 0.30 |
| **B: A + gust forecast (`agent/best_params.json`)** | 72.43 | 73.32 | 72.93 | **73.05** | **+0.67 ± 0.43** |

B is better than A on the development set (+0.32 ± 0.39) and on the test set (+0.37 ± 0.40). Scenarios with B: 49.7,
83.7, 60.8, 61.9 and 67.2 steps. Scenario 5 is worse than with the previous best (59.4).
