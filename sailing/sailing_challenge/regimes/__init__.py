"""
Regimes: the wind generator with its settings. The settings of the public training regime are in train.json.
The validation and test regimes use the same code with hidden settings (Codabench only).
"""
import copy
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PUBLIC_REGIMES = ("train",)


def load_regime(regime="train"):
    """A regime given by its name ("train"), by the path of a JSON file, or as a dict (returned as a copy)."""
    if isinstance(regime, dict):
        return copy.deepcopy(regime)
    if regime in PUBLIC_REGIMES:
        path = os.path.join(HERE, f"{regime}.json")
    elif isinstance(regime, str) and regime.endswith(".json") and os.path.isfile(regime):
        path = regime
    else:
        raise ValueError(f"unknown regime {regime!r}: use one of {PUBLIC_REGIMES} or the path of a JSON file")
    with open(path) as f:
        return json.load(f)


def train_seeds():
    """The fixed public seeds of the training regime, used for the train score."""
    return list(load_regime("train")["fixed_seeds"])
