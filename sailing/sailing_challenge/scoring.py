"""
The score, in a single place (used locally, in the notebooks and on Codabench).

Episode score: 100 * 0.995 ** t if the boat arrives, 0 otherwise, where t is the fractional arrival time:
t = (k - 1) + s if the boat enters the goal circle during step k, after travelling a fraction s of that step.
Score of an agent: mean of the episode scores.
"""
import numpy as np

from .constants import DISCOUNT, GOAL_REWARD


def episode_score(arrival_time):
    return 0.0 if arrival_time is None else float(GOAL_REWARD * DISCOUNT ** arrival_time)


def summarize(results):
    """Aggregate a list of episode results (dicts as returned by evaluation.run_episode)."""
    scores = np.array([r["score"] for r in results], dtype=float)
    times = [r["arrival_time"] for r in results if r["arrival_time"] is not None]
    steps = sum(r["steps"] for r in results)
    return {
        "episodes": len(results),
        "score": float(scores.mean()) if len(scores) else 0.0,
        "score_stderr": float(scores.std(ddof=1) / np.sqrt(len(scores))) if len(scores) > 1 else 0.0,
        "success_rate": len(times) / len(results) if results else 0.0,
        "mean_arrival_time": float(np.mean(times)) if times else None,
        "crashes": int(sum(r["crashes"] for r in results)),
        "invalid_actions": int(sum(r["invalid_actions"] for r in results)),
        "mean_seconds_per_step": sum(r["act_seconds"] for r in results) / steps if steps else 0.0,
        "max_seconds_per_step": max((r.get("max_step_seconds", 0.0) for r in results), default=0.0),
        "errors": sum(r.get("status") == "error" for r in results),
        "timeouts": sum(r.get("status") in ("timeout", "skipped") for r in results),
    }
