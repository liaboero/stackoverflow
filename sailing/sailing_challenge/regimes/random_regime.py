"""
Random regimes: the same wind generator, with settings drawn at random around those of the training regime.

    from sailing_challenge.regimes.random_regime import random_regime
    regime = random_regime(seed)           # a dict, usable as SailingEnv(regime) or run_episode(..., regime=regime)
    env = SailingEnv("random")             # a new random regime (and island layout) at every reset(seed)

Every mechanism keeps its form; what changes is how often and how strongly it acts: direction and variability of
the base wind, period and amplitude of the oscillation, probability, strength and timing of the fronts, number,
size and strength of the gusts, island layout. The draws depend on the seed only.
"""
import copy

import numpy as np

from . import load_regime

LAYOUTS = ("central", "three", "four")


def _range(rng, lo, hi, min_width=0.0):
    """A random sub-range [a, b] of [lo, hi]."""
    a, b = sorted(rng.uniform(lo, hi, 2))
    if b - a < min_width:
        m = 0.5 * (a + b)
        a, b = max(lo, m - 0.5 * min_width), min(hi, m + 0.5 * min_width)
    return [float(a), float(b)]


def random_regime(seed, base="train"):
    rng = np.random.default_rng([int(seed), 7919])
    reg = copy.deepcopy(load_regime(base))
    reg["name"] = f"random-{int(seed)}"
    reg["description"] = "Random regime drawn around the training regime."
    reg.pop("fixed_seeds", None)
    reg["islands"] = LAYOUTS[int(rng.integers(len(LAYOUTS)))]
    w = reg["wind"]
    b = w["base"]
    b["speed"] = float(rng.uniform(8.0, 12.0))
    b["mean_from_deg"] = _range(rng, -140.0, 140.0, 60.0)
    b["spatial_range_deg"] = _range(rng, 10.0, 200.0, 40.0)
    b["mean_from_shape"] = "triangular" if rng.random() < 0.5 else "uniform"
    o = w["oscillation"]
    o["period_steps"] = _range(rng, 25.0, 140.0, 20.0)
    o["amplitude_deg"] = _range(rng, 5.0, 60.0, 10.0)
    o["noise_deg"] = float(rng.uniform(0.5, 4.0))
    f = w["front"]
    f["probability"] = float(rng.uniform(0.2, 0.95))
    f["right_probability"] = float(rng.uniform(0.2, 0.8))
    f["twist_deg"] = _range(rng, 30.0, 80.0, 10.0)
    f["start_y"] = _range(rng, 15.0, 60.0, 10.0)
    f["turning_probability"] = float(rng.uniform(0.2, 0.9))
    f["turning_steps"] = _range(rng, 10.0, 50.0, 10.0)
    g = w["gusts"]
    g["rate_per_step"] = float(rng.uniform(0.02, 0.12))
    g["radius_cells"] = _range(rng, 8.0, 32.0, 5.0)
    g["strength_magnitude"] = _range(rng, 0.2, 0.7, 0.1)
    g["speed_cells_per_step"] = _range(rng, 0.8, 4.0, 0.5)
    g["shift_deg"] = [-float(rng.uniform(5.0, 40.0)), float(rng.uniform(5.0, 40.0))]
    return reg
