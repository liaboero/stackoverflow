"""
Evaluation: exactly the code that runs on Codabench.

    results, summary = evaluate(agent)                 # an agent object, on the 20 fixed seeds of the training regime
    run = evaluate_submission("my_submission/")        # a submission folder or zip, like the Codabench ingestion
    results, summary = score_run(run)                  # replays the recorded actions and scores them, like the
                                                       # Codabench scoring program
    print(format_report(debug_report(results, summary)))

The same rules apply locally and on Codabench:
- the agent is created once (its files are read once); every episode plays a fresh copy of it (copy.deepcopy),
  so nothing is remembered from one episode to the next;
- the agent never sees the seed of the episode;
- an invalid action (None, not an integer, outside 0..8) is replaced by "stay" and counted;
- an exception raised by the agent ends its episode (score 0) and is reported; the next episodes still run;
- time budget: BUDGET_SECONDS for the CODABENCH_EPISODES episodes of a Codabench phase (loading the agent included).
  When it runs out, the current episode stops and the remaining episodes score 0.
  Locally, the budget is scaled to the number of episodes (36 s per episode on average).
- the score is computed by replaying the recorded actions in a fresh environment (scoring step).
"""
import contextlib
import copy
import importlib.util
import math
import os
import re
import signal
import sys
import tempfile
import threading
import time
import traceback
import zipfile

from .constants import ACTION_NAMES, N_ACTIONS, VERSION
from .env import SailingEnv, parse_action
from .regimes import load_regime
from .scoring import episode_score, summarize

BUDGET_SECONDS = 1800.0     # execution_time_limit of a Codabench submission
BUDGET_CLOCK = "wall"       # how the budget is counted: "wall" (real time) or "cpu" (CPU time of the evaluation process)
WALL_LIMIT_FACTOR = 2.0     # with the "cpu" clock, a real-time limit too (this times the budget), for agents that sleep
CLOCKS = {"wall": time.perf_counter, "cpu": time.process_time}
CODABENCH_EPISODES = 50     # episodes of a Codabench phase (validation or test)
MESSAGE_MAX_CHARS = 200     # error messages are truncated to this length in the debug report
AGENT_FILE = "agent.py"
AGENT_CLASS = "MyAgent"
INVALID_ACTION = -1         # how an invalid action is recorded (it is played as "stay")
RUN_FORMAT = 1
_IGNORED = {"__MACOSX", "__pycache__", ".ipynb_checkpoints", ".DS_Store", ".git"}


class TimeBudgetExceeded(BaseException):
    """Raised inside the agent's code when the time budget runs out (a BaseException: `except Exception` misses it)."""

    def __str__(self):
        return "the time budget of the evaluation ran out"


class SubmissionError(Exception):
    """The submission cannot be loaded. `error` describes the exception raised by the agent's code, if any."""

    def __init__(self, message, error=None):
        super().__init__(message if error is None else f"{message}: {format_error(error)}")
        self.reason, self.error = message, error


# ---------------------------------------------------------------------------------------------- errors and limits
def _truncate(text, n=MESSAGE_MAX_CHARS):
    text = str(text)
    return text if len(text) <= n else text[:n - 3] + "..."


def describe_error(exc, source_dir=None):
    """Type, location in the submission's files and message (truncated) of an exception raised by an agent."""
    try:
        message = str(exc)
    except Exception:  # a broken __str__ in the agent's own exception
        message = ""
    location = None
    if source_dir is not None:
        root = os.path.realpath(source_dir)
        frames = [(f.filename, f.lineno, f.name) for f in traceback.extract_tb(exc.__traceback__)]
        if isinstance(exc, SyntaxError) and exc.filename:
            frames.append((exc.filename, exc.lineno, "<module>"))
        for filename, lineno, name in reversed(frames):
            path = os.path.realpath(filename)
            if path.startswith(root + os.sep):
                location = f"{os.path.relpath(path, root)}, line {lineno}, in {name}"
                break
    return {"type": _truncate(type(exc).__name__, 60), "message": _truncate(message),
            "location": None if location is None else _truncate(location)}


def _mask_digits(text):
    """Digits replaced by '#': an error message must not carry numbers out of a hidden race (question 26)."""
    return re.sub(r"[0-9]", "#", str(text))


def _mask_location(location):
    """'file, line 8, in name': the line number is kept (it only depends on the code), the rest has its digits masked."""
    m = re.fullmatch(r"(.*), line ([0-9]+), in (.*)", str(location), flags=re.S)
    if m is None:
        return _mask_digits(location)
    return f"{_mask_digits(m.group(1))}, line {m.group(2)}, in {_mask_digits(m.group(3))}"


def sanitize_error(error):
    """
    The fields of an error description, with their types and lengths enforced (the scoring step trusts nothing),
    and the digits of the type, message and location masked (the line number of the location is kept).
    """
    if not isinstance(error, dict):
        return None
    location = error.get("location")
    return {"type": _mask_digits(_truncate(error.get("type", ""), 60)), "message": _mask_digits(_truncate(error.get("message", ""))),
            "location": None if location is None else _mask_location(_truncate(location))}


def format_error(error):
    where = f" ({error['location']})" if error.get("location") else ""
    message = f": {error['message']}" if error.get("message") else ""
    return f"{error['type']}{where}{message}"


def describe_invalid_action(action):
    """What was wrong with an invalid action (its type only: the debug report never shows an arbitrary value)."""
    if action is None:
        return "None"
    if isinstance(action, (bool,)) or type(action).__name__ == "bool_":
        return "a boolean"
    a, valid = parse_action(action)
    if not valid:
        try:
            if isinstance(action, int) or (hasattr(action, "dtype") and getattr(action, "ndim", 1) == 0
                                           and action.dtype.kind in "iu"):
                return f"the integer {int(action)}, outside 0..{N_ACTIONS - 1}"
        except Exception:
            pass
    shape = getattr(action, "shape", None)
    kind = _truncate(type(action).__name__, 40)
    return f"a {kind} of shape {tuple(shape)}" if shape is not None else f"a {kind}"


class _Deadline:
    """
    Time budget. clock="wall": real time (SIGALRM). clock="cpu": CPU time of the process, all threads (SIGPROF), so that
    an agent is not penalised when the machine is shared; plus a real-time limit (wall_limit seconds, by default
    WALL_LIMIT_FACTOR times the budget), so that an agent that sleeps or waits still stops. With the signals (Linux,
    macOS, main thread), an agent stuck in act() is interrupted, and the report tells where it was.
    """

    def __init__(self, seconds, clock="wall", wall_limit=None):
        if clock not in CLOCKS:
            raise ValueError(f"clock must be one of {sorted(CLOCKS)}, not {clock!r}")
        self.clock_name, self.clock = clock, CLOCKS[clock]
        self.end = None if seconds is None else self.clock() + max(float(seconds), 0.0)
        if clock == "cpu" and seconds is not None:
            wall = WALL_LIMIT_FACTOR * max(float(seconds), 0.0) if wall_limit is None else max(float(wall_limit), 0.0)
            self.wall_end = time.perf_counter() + wall
        else:
            self.wall_end = None
        self._timers = []

    def _arm(self, sig, timer, seconds):
        previous = signal.signal(sig, self._ring)
        signal.setitimer(timer, max(seconds, 1e-3))
        self._timers.append((sig, timer, previous))

    def __enter__(self):
        if (self.end is not None and hasattr(signal, "setitimer")
                and threading.current_thread() is threading.main_thread()):
            if self.clock_name == "wall":
                self._arm(signal.SIGALRM, signal.ITIMER_REAL, self.end - time.perf_counter())
            else:
                self._arm(signal.SIGPROF, signal.ITIMER_PROF, self.end - time.process_time())
                self._arm(signal.SIGALRM, signal.ITIMER_REAL, self.wall_end - time.perf_counter())
        return self

    def __exit__(self, *exc):
        while self._timers:
            sig, timer, previous = self._timers.pop()
            signal.setitimer(timer, 0)
            signal.signal(sig, previous)
        return False

    @staticmethod
    def _ring(signum, frame):
        raise TimeBudgetExceeded()

    def remaining(self):
        """(seconds left on the budget clock, real seconds left or None): for a nested deadline."""
        if self.end is None:
            return None, None
        wall = None if self.wall_end is None else self.wall_end - time.perf_counter()
        return self.end - self.clock(), wall

    def expired(self):
        if self.end is None:
            return False
        return self.clock() >= self.end or (self.wall_end is not None and time.perf_counter() >= self.wall_end)

    # The signals can ring a few milliseconds late (about 10 ms on macOS, timer coalescing; SIGPROF is sampled at the
    # scheduler tick). While they are armed, the check at the start of each step leaves them this margin: the signal,
    # which tells where the agent was, rings first.
    ALARM_GRACE = 0.05

    def check(self):
        if self.end is None:
            return
        grace = self.ALARM_GRACE if self._timers else 0.0
        if self.clock() >= self.end + grace or (self.wall_end is not None and time.perf_counter() >= self.wall_end + grace):
            raise TimeBudgetExceeded()


# ---------------------------------------------------------------------------------------------- episodes
def _failed_episode(seed, status, error=None):
    return {"seed": seed, "arrived": False, "arrival_time": None, "score": 0.0, "steps": 0, "crashes": 0,
            "invalid_actions": 0, "act_seconds": 0.0, "max_step_seconds": 0.0, "status": status, "error": error,
            "error_step": None, "first_invalid_action": None, "actions": [], "trajectory": None}


def run_episode(agent, seed, regime="train", static_wind=False, physics=None, record=False, env=None,
                catch_errors=False, deadline=None, source_dir=None):
    """
    Run one episode and return its result (dict). The agent never sees the seed.

    catch_errors: if True, an exception raised by the agent ends the episode (score 0) and is described in
    result["error"]; if False (default), it propagates. `deadline` and `source_dir` are used by evaluate().
    """
    env = env if env is not None else SailingEnv(regime, static_wind=static_wind, physics=physics)
    clock = deadline.clock if deadline is not None else time.perf_counter   # act() is timed on the budget's clock
    obs, info = env.reset(seed=seed)
    trajectory = [obs["position"].tolist()] if record else None
    actions, act_seconds, max_step_seconds = [], 0.0, 0.0
    status, error, first_invalid = "ok", None, None
    try:
        agent.reset(info["env_info"])
        terminated = truncated = False
        while not (terminated or truncated):
            if deadline is not None:
                deadline.check()
            t0 = clock()
            action = agent.act(obs)
            dt = clock() - t0
            act_seconds += dt
            max_step_seconds = max(max_step_seconds, dt)
            a, valid = parse_action(action)
            if not valid and first_invalid is None:
                first_invalid = describe_invalid_action(action)
            actions.append(a if valid else INVALID_ACTION)  # before the step: the replay plays it in any case
            obs, _, terminated, truncated, info = env.step(action)
            if record:
                trajectory.append(obs["position"].tolist())
    except TimeBudgetExceeded as exc:
        status, error = "timeout", describe_error(exc, source_dir)
    except (Exception, SystemExit) as exc:
        if not catch_errors:
            raise
        status, error = "error", describe_error(exc, source_dir)
    return {
        "seed": seed,
        "arrived": env.arrival_time is not None,
        "arrival_time": env.arrival_time,
        "score": episode_score(env.arrival_time),
        "steps": env.t,
        "crashes": env.n_crashes,
        "invalid_actions": env.n_invalid_actions,
        "act_seconds": act_seconds,
        "max_step_seconds": max_step_seconds,
        "status": status,
        "error": error,
        "error_step": env.t if error is not None else None,
        "first_invalid_action": first_invalid,
        "actions": actions,
        "trajectory": trajectory,
    }


def evaluate(agent, seeds=None, regime="train", static_wind=False, physics=None, record=False, *, fresh_copy=True,
             catch_errors=True, time_budget=None, source_dir=None, progress=None, clock="wall", wall_limit=None):
    """
    Run the agent on several seeds, as on Codabench, and return (results, summary).

    seeds: the fixed seeds of the regime by default (for "train": the 20 public seeds, see train_seeds()).
    fresh_copy: every episode plays copy.deepcopy(agent); the agent itself is never played.
    catch_errors: an exception of the agent ends its episode (score 0) instead of stopping the evaluation.
    time_budget: seconds for the whole evaluation (None: no limit). The remaining episodes score 0 when it runs out.
    progress: optional function called with (index, result) after every episode.
    clock: how the budget is counted, "wall" (real time) or "cpu" (CPU time, with a real-time limit wall_limit).
    """
    if seeds is None:
        seeds = load_regime(regime)["fixed_seeds"]
    seeds = [int(s) for s in seeds]
    env = SailingEnv(regime, static_wind=static_wind, physics=physics)
    results = []
    start, wall_start = CLOCKS[clock](), time.perf_counter()
    deadline = _Deadline(time_budget, clock, wall_limit)
    try:
        with deadline:
            for seed in seeds:
                if deadline.expired():
                    break
                try:
                    player = copy.deepcopy(agent) if fresh_copy else agent
                except TimeBudgetExceeded:
                    break
                except (Exception, SystemExit) as exc:
                    if not catch_errors:
                        raise
                    error = describe_error(exc, source_dir)
                    error["message"] = _truncate(f"copy.deepcopy(agent) failed: {error['message']}")
                    results.append(_failed_episode(seed, "error", error))
                    continue
                results.append(run_episode(player, seed, env=env, record=record, catch_errors=catch_errors,
                                           deadline=deadline if time_budget is not None else None,
                                           source_dir=source_dir))
                results[-1]["episode"] = len(results)
                if progress is not None:
                    progress(len(results) - 1, results[-1])
    except TimeBudgetExceeded:  # rang between two episodes
        pass
    for seed in seeds[len(results):]:
        results.append(_failed_episode(seed, "skipped"))
    for i, r in enumerate(results):
        r["episode"] = i + 1
    summary = summarize(results)
    summary["total_seconds"] = CLOCKS[clock]() - start
    summary["wall_seconds"] = time.perf_counter() - wall_start
    summary["budget_seconds"] = time_budget
    summary["clock"] = clock
    return results, summary


# ---------------------------------------------------------------------------------------------- submissions
_loaded_roots = []


def find_agent_root(folder):
    """The folder that holds agent.py: the folder itself, or its only subfolder (a zip of the folder)."""
    if os.path.isfile(os.path.join(folder, AGENT_FILE)):
        return folder
    entries = [e for e in os.listdir(folder) if e not in _IGNORED and not e.startswith(".")]
    if len(entries) == 1 and os.path.isfile(os.path.join(folder, entries[0], AGENT_FILE)):
        return os.path.join(folder, entries[0])
    raise SubmissionError(f"{AGENT_FILE} not found at the root of the submission (found: {sorted(entries)[:10]})")


def _forget_previous_submissions():
    """Unload the modules of a submission loaded earlier in this process (local use: several submissions in turn)."""
    for root in _loaded_roots:
        while root in sys.path:
            sys.path.remove(root)
        for name, module in list(sys.modules.items()):
            path = getattr(module, "__file__", None)
            if path and os.path.realpath(path).startswith(root + os.sep):
                del sys.modules[name]
    _loaded_roots.clear()
    sys.modules.pop("agent", None)


def load_agent(folder):
    """
    Import agent.py from a submission folder and create MyAgent() once. Returns (agent, root folder).
    The folder is put first on sys.path, so agent.py can import its own helper modules.
    Raises SubmissionError with a clear message if anything fails.
    """
    import sailing_challenge  # noqa: F401  (already imported: a folder named sailing_challenge cannot replace it)
    root = os.path.realpath(find_agent_root(folder))
    _forget_previous_submissions()
    sys.path.insert(0, root)
    _loaded_roots.append(root)
    spec = importlib.util.spec_from_file_location("agent", os.path.join(root, AGENT_FILE))
    module = importlib.util.module_from_spec(spec)
    sys.modules["agent"] = module
    try:
        spec.loader.exec_module(module)
    except (Exception, SystemExit) as exc:
        raise SubmissionError(f"importing {AGENT_FILE} failed", describe_error(exc, root)) from None
    cls = getattr(module, AGENT_CLASS, None)
    if not isinstance(cls, type):
        raise SubmissionError(f"{AGENT_FILE} must define a class {AGENT_CLASS} (a subclass of BaseAgent)")
    try:
        agent = cls()
    except (Exception, SystemExit) as exc:
        raise SubmissionError(f"{AGENT_CLASS}() failed", describe_error(exc, root)) from None
    for method in ("reset", "act"):
        if not callable(getattr(agent, method, None)):
            raise SubmissionError(f"{AGENT_CLASS} has no method {method}(): derive it from BaseAgent")
    return agent, root


@contextlib.contextmanager
def submission_folder(submission):
    """A submission given as a folder or a zip file: yields a folder (the zip is extracted to a temporary folder).

    The folder is absolute: the agent then runs from another working directory."""
    submission = os.path.abspath(os.fspath(submission))
    if os.path.isdir(submission):
        yield submission
        return
    if not zipfile.is_zipfile(submission):
        raise SubmissionError(f"{submission!r} is neither a folder nor a zip file")
    with tempfile.TemporaryDirectory(prefix="submission_") as tmp:
        with zipfile.ZipFile(submission) as z:
            z.extractall(tmp)
        yield tmp


@contextlib.contextmanager
def _working_directory():
    """Run the agent from an empty temporary folder, as on Codabench: files must be found relative to __file__."""
    previous = os.getcwd()
    with tempfile.TemporaryDirectory(prefix="cwd_") as tmp:
        os.chdir(tmp)
        try:
            yield tmp
        finally:
            os.chdir(previous)


def default_budget(n_episodes):
    """The local budget: BUDGET_SECONDS scaled to the number of episodes."""
    return BUDGET_SECONDS * n_episodes / CODABENCH_EPISODES


def evaluate_submission(submission, regime="train", seeds=None, time_budget="auto", progress=None, clock=None,
                        wall_limit=None):
    """
    Evaluate a submission (folder or zip) like the Codabench ingestion step, and return the run (a JSON-ready dict:
    for every episode the recorded actions, the status, the error if any and the timing). Score it with score_run().
    time_budget: "auto" = BUDGET_SECONDS scaled to the number of seeds; None = no limit.
    clock: "wall" or "cpu" (BUDGET_CLOCK by default); wall_limit: the real-time limit with the "cpu" clock (seconds).
    """
    clock = BUDGET_CLOCK if clock is None else clock
    regime = load_regime(regime)
    seeds = [int(s) for s in (regime["fixed_seeds"] if seeds is None else seeds)]
    if time_budget == "auto":
        time_budget = default_budget(len(seeds))
    run = {"format": RUN_FORMAT, "version": VERSION, "regime": regime.get("name"), "seeds": seeds,
           "budget_seconds": time_budget, "clock": clock, "load_seconds": None, "load_error": None,
           "total_seconds": None, "wall_seconds": None, "episodes": []}
    start, wall_start = CLOCKS[clock](), time.perf_counter()
    deadline = _Deadline(time_budget, clock, wall_limit)
    with submission_folder(submission) as folder, _working_directory():
        try:
            with deadline:
                agent, root = load_agent(folder)
        except SubmissionError as exc:
            run["load_error"] = {"type": "SubmissionError", "message": _truncate(exc.reason), "location": None} \
                if exc.error is None else dict(exc.error, message=_truncate(f"{exc.reason}: {exc.error['message']}"))
        except TimeBudgetExceeded as exc:
            run["load_error"] = describe_error(exc)
        run["load_seconds"] = CLOCKS[clock]() - start
        if run["load_error"] is None:
            remaining, wall_remaining = deadline.remaining()
            results, _ = evaluate(agent, seeds, regime, time_budget=remaining, source_dir=root, progress=progress,
                                  clock=clock, wall_limit=wall_remaining)
            keep = ("episode", "seed", "status", "error", "error_step", "first_invalid_action", "actions",
                    "steps", "act_seconds", "max_step_seconds")
            run["episodes"] = [{k: r[k] for k in keep} for r in results]
    run["total_seconds"] = CLOCKS[clock]() - start
    run["wall_seconds"] = time.perf_counter() - wall_start
    return run


# ---------------------------------------------------------------------------------------------- scoring
def _number(x, default=0.0):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return default
    return x if math.isfinite(x) and x >= 0 else default


def replay_episode(seed, actions, regime="train", env=None):
    """Play recorded actions (ints, INVALID_ACTION for an invalid one) in a fresh episode; return its outcome."""
    env = env if env is not None else SailingEnv(regime)
    env.reset(seed=seed)
    done = False
    for a in actions if isinstance(actions, list) else []:
        if done:
            break
        a = a if isinstance(a, int) and not isinstance(a, bool) and 0 <= a < N_ACTIONS else INVALID_ACTION
        *_, terminated, truncated, _ = env.step(a)
        done = terminated or truncated
    return {"seed": seed, "arrived": env.arrival_time is not None, "arrival_time": env.arrival_time,
            "score": episode_score(env.arrival_time), "steps": env.t, "crashes": env.n_crashes,
            "invalid_actions": env.n_invalid_actions}


def score_run(run, regime="train", seeds=None):
    """
    Score a run (the output of evaluate_submission) like the Codabench scoring step, and return (results, summary).
    The outcome of every episode is recomputed by replaying its actions; the seeds come from the regime, not from
    the run, and a missing episode scores 0. Only the timing and the error descriptions are taken from the run.
    """
    regime = load_regime(regime)
    seeds = [int(s) for s in (regime["fixed_seeds"] if seeds is None else seeds)]
    run = run if isinstance(run, dict) else {}
    episodes = {}
    for ep in run.get("episodes", []) if isinstance(run.get("episodes"), list) else []:
        if isinstance(ep, dict) and isinstance(ep.get("seed"), int):
            episodes.setdefault(ep["seed"], ep)
    load_error = sanitize_error(run.get("load_error"))
    env = SailingEnv(regime)
    results = []
    for i, seed in enumerate(seeds):
        ep = episodes.get(seed)
        if ep is None:
            r = _failed_episode(seed, "not_run" if load_error is not None else "missing")
        else:
            r = replay_episode(seed, ep.get("actions"), env=env)
            status = ep.get("status")
            r.update({
                "status": status if status in ("ok", "error", "timeout", "skipped") else "ok",
                "error": sanitize_error(ep.get("error")),
                "error_step": int(ep["error_step"]) if isinstance(ep.get("error_step"), int) else None,
                "first_invalid_action": None if ep.get("first_invalid_action") is None
                else _truncate(ep["first_invalid_action"], 80),
                "act_seconds": _number(ep.get("act_seconds")),
                "max_step_seconds": _number(ep.get("max_step_seconds")),
            })
        r["episode"] = i + 1
        results.append(r)
    summary = summarize(results)
    clock = run.get("clock") if run.get("clock") in CLOCKS else "wall"   # runs written before the clock option: wall
    summary.update({"total_seconds": _number(run.get("total_seconds"), None),
                    "wall_seconds": _number(run.get("wall_seconds"), None), "clock": clock,
                    "budget_seconds": _number(run.get("budget_seconds"), None),
                    "load_seconds": _number(run.get("load_seconds"), None), "load_error": load_error})
    return results, summary


# ---------------------------------------------------------------------------------------------- debug report
def debug_report(results, summary):
    """
    The debug report shown to the participants (also for the hidden regimes): errors with their type, location in
    the submission's files and message (truncated), episode concerned, time budget, invalid actions, time per step.
    Never a position, a trajectory or a wind.
    """
    problems = []
    for r in results:
        if r.get("status", "ok") not in ("error", "timeout") and not r.get("invalid_actions"):
            continue  # the episodes not played are only counted
        steps = r.get("steps") or 0
        problems.append({
            "episode": r["episode"],
            "status": r.get("status", "ok"),
            "error": sanitize_error(r.get("error")),
            "error_step": r.get("error_step"),
            "invalid_actions": int(r.get("invalid_actions", 0)),
            "first_invalid_action": r.get("first_invalid_action"),
            "mean_ms_per_step": 1000 * _number(r.get("act_seconds")) / steps if steps else 0.0,
            "max_ms_per_step": 1000 * _number(r.get("max_step_seconds")),
        })
    statuses = [r.get("status", "ok") for r in results]
    completed_scores = [_number(r.get("score")) for r in results if r.get("status", "ok") == "ok"]
    return {
        "score": summary.get("score"),
        "completed_score": sum(completed_scores) / len(completed_scores) if completed_scores else None,
        "episodes": len(results),
        "completed": statuses.count("ok"),
        "errors": statuses.count("error"),
        "timeouts": statuses.count("timeout"),
        "skipped": statuses.count("skipped"),
        "not_run": statuses.count("not_run") + statuses.count("missing"),
        "first_not_played": min((r["episode"] for r in results if r.get("status") in ("skipped", "not_run", "missing")),
                                default=None),
        "invalid_actions": int(summary.get("invalid_actions", 0)),
        "load_error": sanitize_error(summary.get("load_error")),
        "load_seconds": summary.get("load_seconds"),
        "total_seconds": summary.get("total_seconds"),
        "wall_seconds": summary.get("wall_seconds"),
        "clock": summary.get("clock", "wall"),
        "budget_seconds": summary.get("budget_seconds"),
        "mean_ms_per_step": 1000 * summary.get("mean_seconds_per_step", 0.0),
        "max_ms_per_step": 1000 * summary.get("max_seconds_per_step", 0.0),
        "problems": problems,
    }


_STATUS_TEXT = {"error": "stopped by an exception", "timeout": "stopped: time budget exceeded",
                "skipped": "not played: time budget exceeded", "ok": "completed"}


def report_warnings(report):
    """The warnings shown first, in capitals, when the score is not the level of the agent: races not played or
    stopped (they score 0). A list of (headline, explanation)."""
    out = []
    if report["load_error"] is not None:
        out.append(("WARNING: YOUR SUBMISSION COULD NOT BE LOADED. NO RACE WAS PLAYED: YOUR SCORE IS 0.",
                    "The error is below. Check your zip locally: python -m sailing_challenge.submission check submission.zip"))
        return out
    n, done = report["episodes"], report["completed"]
    over = report["timeouts"] + report["skipped"]
    level = ""
    if report.get("completed_score") is not None and report.get("score") is not None and done < n:
        level = (f" Your score ({report['score']:.2f}) is therefore NOT the level of your agent: on the {done} race(s) "
                 f"it completed, it scores {report['completed_score']:.2f} on average.")
    if over:
        budget = report.get("budget_seconds")
        unit = " of CPU time" if report.get("clock") == "cpu" else ""
        within = f" within the time budget ({budget:.0f} s{unit})" if budget is not None else " within the time budget"
        out.append(("WARNING: YOUR AGENT DID NOT FINISH ALL THE RACES WITHIN THE TIME BUDGET.",
                    f"Only {done} of the {n} races were completed{within}: the {over} other(s) score 0.{level} "
                    "Make act() faster (see the time per step below), and check the projected time locally: "
                    "python -m sailing_challenge.submission check submission.zip"))
        level = ""
    if report["errors"]:
        out.append((f"WARNING: {report['errors']} RACE(S) ENDED WITH AN EXCEPTION OF YOUR AGENT, AND SCORE 0.",
                    f"The error and where it happened are below.{level}"))
        level = ""
    if report["not_run"]:
        out.append((f"WARNING: THE EVALUATION WAS INTERRUPTED: {report['not_run']} RACE(S) HAVE NO RESULT, AND SCORE 0.",
                    f"Details below.{level}"))
    return out


def format_report(report, warnings=True):
    """The debug report as text: the warnings first, if any (warnings=False: without them, shown apart)."""
    lines = []
    for headline, explanation in (report_warnings(report) if warnings else []):
        lines += [headline, explanation, ""]
    if report["load_error"] is not None:
        lines.append(f"The submission could not be loaded: {format_error(report['load_error'])}")
        lines.append("No episode was played (score 0).")
        return "\n".join(lines)
    total, budget = report.get("total_seconds"), report.get("budget_seconds")
    cpu = report.get("clock") == "cpu"
    timing = f"{total:.1f} s" if total is not None else "unknown time"
    if budget is not None:
        timing += f" (budget {budget:.0f} s{' of CPU time' if cpu else ''})"
    if cpu and report.get("wall_seconds") is not None:
        timing += f", {report['wall_seconds']:.1f} s of real time"
    lines.append(f"Episodes: {report['episodes']}, completed: {report['completed']}, errors: {report['errors']}, "
                 f"time budget exceeded: {report['timeouts'] + report['skipped']}, "
                 f"invalid actions: {report['invalid_actions']}.")
    lines.append(f"{'Total CPU time' if cpu else 'Total time'}: {timing}. "
                 f"{'CPU time' if cpu else 'Time'} per step: {report['mean_ms_per_step']:.2f} ms on average, "
                 f"{report['max_ms_per_step']:.1f} ms at most.")
    if report["not_run"]:
        lines.append(f"{report['not_run']} episode(s) have no result, from episode {report['first_not_played']}: "
                     "the evaluation was interrupted (time limit of the platform, memory, or the process exited). "
                     "They score 0.")
    if report["timeouts"] or report["skipped"]:
        lines.append(f"The time budget ran out: {report['skipped']} episode(s) were not played, from episode "
                     f"{report['first_not_played']}, and score 0. Make act() faster.")
    groups = {}  # identical problems in several episodes are shown once
    for p in report["problems"]:
        key = (p["status"], None if p["error"] is None else format_error(p["error"]), p["first_invalid_action"])
        groups.setdefault(key, []).append(p)
    for (status, error, first_invalid), problems in groups.items():
        p = problems[0]
        episodes = [q["episode"] for q in problems]
        head = (f"Episode {p['episode']}" if len(problems) == 1
                else f"Episodes {_ranges(episodes)} ({len(problems)} episodes, details of episode {p['episode']})")
        head += f": {_STATUS_TEXT.get(status, status)}"
        if p["error_step"] is not None and status in ("error", "timeout"):
            head += f" at step {p['error_step']}"
        if p["error"] is not None and status == "error":
            head += f": {error}"
        elif p["error"] is not None and status == "timeout" and p["error"].get("location"):
            head += f" (the agent was in {p['error']['location']})"
        if p["invalid_actions"]:
            head += (f"; {p['invalid_actions']} invalid action(s) played as STAY, the first one was "
                     f"{first_invalid} (expected an int in 0..{N_ACTIONS - 1}: "
                     f"{', '.join(f'{i} = {n}' for i, n in enumerate(ACTION_NAMES))})")
        head += f"; {p['mean_ms_per_step']:.2f} ms per step on average, {p['max_ms_per_step']:.1f} ms at most"
        lines.append(head + ".")
    if not report["problems"]:
        lines.append("No error, no invalid action.")
    return "\n".join(lines)


def _ranges(numbers):
    """[1, 2, 3, 7] -> "1-3, 7"."""
    parts, start = [], None
    for i, n in enumerate(numbers):
        if start is None:
            start = n
        if i + 1 == len(numbers) or numbers[i + 1] != n + 1:
            parts.append(str(start) if start == n else f"{start}-{n}")
            start = None
    return ", ".join(parts)
