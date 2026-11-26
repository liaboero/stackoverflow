# Sailing challenge: agents for wind maps never seen before

Test set: 75 episodes never used in any training. It holds the 5 game scenarios, 30 new seeds of the training
regime and 40 random regimes (wind settings and island layout drawn by `regimes/random_regime.py`). The table gives
the mean score, 100 · 0.995^steps.

| Agent | Scenarios | Training regime | Random regimes | All |
|---|---|---|---|---|
| ForecastPlanner, trained parameters (`agent/`) | 72.65 | 73.48 | 70.96 | **72.08** |
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
