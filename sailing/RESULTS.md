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

## Retraining with the new forecasts (`agent/best_params_retrained.json`)

8 parameters were retrained jointly, starting from B: 10 generations on 64 episodes, training score 74.09 → 74.75.

| Version | Scenarios | Training regime | Random regimes | All (test) | Development set |
|---|---|---|---|---|---|
| B (`best_params.json`) | 72.43 | 73.32 | 72.93 | 73.05 | 72.53 |
| Retrained | **73.61** | **73.78** | 72.40 | 73.03 | 72.26 |

Paired difference, retrained - B: -0.02 ± 0.34 on the test set, -0.26 ± 0.28 on the development set. In general the two
are equivalent: the retraining brings no measurable gain. On the 5 game scenarios (only 5 races), the retrained
version does 51.9, 78.6, 61.6, 62.6 and 51.9 steps. That beats the game's Forecaster (72.76 on average), and beats the
game's Oracle on scenarios 1 and 2.

## Decisions on 400-episode sets (`agent/validate.py`, `agent/train_robust.py`)

Validation set: 400 new episodes, 300 on random regimes and 100 on the training regime. A paired difference is
measured to about ±0.15 point.

| Version | Validation | Difference vs B |
|---|---|---|
| B (`best_params.json`) | 72.73 | — |
| Retrained | 72.78 | +0.05 ± 0.16 |
| Previous best, without the new forecasts | 72.17 | **-0.56 ± 0.19** |

The gust and oscillation forecasts are a real gain. The retraining is not.

**Training that resists over-fitting.** Every generation sails fresh episodes, 3/4 of them on random regimes. A
candidate must beat the incumbent on a first batch, then again on a second fresh batch. In 10 generations, 5
changes were accepted (generations 2, 3, 5, 6 and 9). The best intermediate version on the validation set was
generation 3, at +0.24 ± 0.11. On a second set of 400 new episodes it does **-0.03 ± 0.11**: its lead was the luck of
being chosen among 5 versions on the same episodes.

**Conclusion:** B stays the best model for unknown maps (72.7 to 72.8 on 800 new episodes). Its parameters are on a
plateau: no retraining improves them measurably. Further gains need new forecasting ideas, not more tuning.

## Scenario 5: false start, cause and fix (`gust_max_tau`)

With B, the boat went west for 6 steps at the start of scenario 5, then came back east (67.2 steps; 57.3 without the
gust forecast). The gust forecast was right in itself: the bottom of the map on the east side loses about 10% of its
wind 30 steps later. But the time-to-goal map uses one forecast field, about 30-40 steps ahead, for the whole route.
The boat passes the bottom of the map in the first 5-15 steps. A moving gust forecast 30 steps ahead was therefore
applied where the boat would be much earlier: a false lull on the good side, worth about 10 steps on the map.

Fix: the gust forecast stops at 10 steps (`gust_max_tau`). Scenario 5: 53.6 steps. On the 400 validation episodes:
+0.04 ± 0.10 (a cap of 20: +0.11 ± 0.09). The fix is neutral on average and corrects a real modelling flaw, so it is
now in `best_params.json` (B without the cap: `best_params_B_nocap.json`). With it, the 5 scenarios take 51.95, 84.51,
60.83, 61.92 and 53.60 steps.

A second flaw was found: the forecast horizon of the map flips from step to step (28, 40, 28 steps...), because it
is set from the remaining time read on the previous map. Making it stable (a map of the present wind every step) is
worse on scenario 5 (81.8 steps). It is left off (`lag_from`).

## False routes on the training episodes, and what they say about the learning (`agent/false_routes.py`)

Reference: the same planner, given the true future wind ("all true"). Results per episode in
`agent/false_routes.json`.

| | Training episodes (64) | Unseen episodes (48) |
|---|---|---|
| Mean loss vs "all true" | **-0.04 steps** | **+2.84 steps** |
| Episodes losing > 5 steps | 9 (14%) | 10 (21%) |
| Share of the lost time from those | 63% | 84% |
| Other side of the central island than "all true" | 13 | 15 |

The largest false routes on the training episodes, and the variant that fixes each one:

| Episode | Model | All true | Cause found |
|---|---|---|---|
| train 1105 (no front) | 69.0 | 45.1 | the agent believes in a turning front that does not exist: front forecast off → 45.2 |
| train 1000 (no front) | 77.7 | 61.7 | oscillation forecast: true oscillation → 61.7 (oscillation forecast off: 95.5) |
| train 1001 (no front) | 93.3 | 79.3 | mixed: oscillation forecast off → 65.5, true gusts → 71.1 |
| random 10005 (no front) | 74.4 | 66.6 | gusts: true gusts → 68.1 |
| train 1112 (turning front) | 73.0 | 65.8 | gusts (true gusts → 64.7) or front (front forecast off → 67.8) |

Conclusions:

1. **A false route is a wrong side.** The boat goes around the central island on the wrong side. The choice is made
   early, is binary, and cannot be undone. A few episodes carry most of the lost time.
2. **The causes vary:** a false turning front (the agent sees one in a base wind pattern), an oscillation forecast
   that is off, a gust forecast. No single parameter fixes them all. Moving one parameter flips some side choices one
   way and others the other way (the gust cap: 16 races better, 10 worse). That is why the tuning reached a plateau.
3. **The training over-fits these coin flips.** On its own training episodes, the model is as good as the planner
   with the true wind (-0.04 steps). It even beats it by more than 5 steps on 7 episodes (115 steps in total), where
   its forecast errors happen to help. On unseen episodes the same gap is +2.84 steps. The parameters have learned
   which way to flip the training episodes' coins, not how to forecast better. With 32-64 episodes, the mean score
   is driven by about ten binary decisions: far too few to learn from.
4. **The planner itself also makes false routes.** With the true future wind, it takes 107 steps on train 1003, where
   the model takes 63.6. Its single-horizon map applies the wind of one future time to the whole route. The flaw
   found on the gusts (scenario 5) is general.

What this means for the learning: more episodes per decision (hundreds), and a signal closer to the decision
itself (was the side right?) rather than only the final score. Above all, structural changes rather than tuning:
choose the side with the uncertainty of the forecast in mind (expected time over several plausible futures, or wait
before committing), detect the front with a probability rather than a threshold, and replace the single-horizon
map.
