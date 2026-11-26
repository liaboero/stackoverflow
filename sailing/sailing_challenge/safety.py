"""
Safety helpers. They are exact: the next step only depends on the observed position, velocity and local wind.

    is_safe(obs, action, env_info)  -> True if the action does not crash the boat on the next step
    safe_actions(obs, env_info)     -> the actions that do not crash it
"""
from .constants import N_ACTIONS
from .dynamics import transition
from .physics import Physics


def next_step(obs, action, env_info):
    """What the environment will do on the next step with this action (a dynamics.Transition)."""
    return transition(obs["position"], obs["velocity"], obs["local_wind"], action, env_info["island_mask"],
                      Physics(**env_info["physics"]), env_info["goal"], env_info["goal_radius"],
                      env_info["grid_size"])


def is_safe(obs, action, env_info):
    return next_step(obs, action, env_info).event != "crash"


def safe_actions(obs, env_info):
    return [a for a in range(N_ACTIONS) if is_safe(obs, a, env_info)]
