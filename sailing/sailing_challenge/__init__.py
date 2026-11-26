"""Sailing challenge (v3): environment, boat physics, safety helpers, distance around the islands, score and local evaluation."""
from .agents import BaseAgent
from .constants import (ACTION_NAMES, ACTION_VECTORS, DISCOUNT, GOAL, GOAL_RADIUS, GRID_SIZE, MAX_STEPS, START,
                        STAY, VERSION)
from .distance import distance_to_goal
from .env import SailingEnv
from .evaluation import evaluate, run_episode
from .evaluation import (BUDGET_SECONDS, SubmissionError, debug_report, evaluate_submission, format_report,
                         load_agent, replay_episode, score_run)
from .physics import Physics
from .regimes import load_regime, train_seeds
from .replay import build_replay, export_replays, load_replay, record_replay, save_replay
from .safety import is_safe, next_step, safe_actions
from .scoring import episode_score

__version__ = VERSION
__all__ = ["ACTION_NAMES", "ACTION_VECTORS", "BaseAgent", "DISCOUNT", "GOAL", "GOAL_RADIUS", "GRID_SIZE",
           "MAX_STEPS", "Physics", "START", "STAY", "SailingEnv", "VERSION", "distance_to_goal", "episode_score", "evaluate",
           "is_safe", "load_regime", "next_step", "run_episode", "safe_actions", "train_seeds"]
