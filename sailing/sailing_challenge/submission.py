"""
Prepare, check and evaluate a submission, exactly as Codabench will run it.

A submission is a zip with agent.py at its root, which defines a class MyAgent(BaseAgent), plus any files it needs
(weights in .npz, scikit-learn models in a pickle, helper modules). Load them with a path relative to __file__:
the evaluation does not run from your folder. See submission_template/ in the repository.

    python -m sailing_challenge.submission zip my_agent/ -o submission.zip   # build the zip, then check it
    python -m sailing_challenge.submission check submission.zip              # check a zip (or a folder)
    python -m sailing_challenge.submission evaluate submission.zip           # full local evaluation (train seeds)

The same functions are available from Python: make_submission_zip, check_submission, evaluate_submission.
"""
import argparse
import copy
import json
import os
import sys
import time
import zipfile

from .agents import BaseAgent
from .constants import ACTION_NAMES, N_ACTIONS
from .env import SailingEnv, parse_action
from .evaluation import (_IGNORED, AGENT_FILE, AGENT_CLASS, BUDGET_SECONDS, CODABENCH_EPISODES, SubmissionError,
                         _working_directory, debug_report, describe_error, describe_invalid_action,
                         evaluate_submission, format_error, format_report, load_agent, score_run,
                         submission_folder)
from .regimes import train_seeds

MAX_ZIP_MB = 100          # a larger zip is reported as a warning
TYPICAL_STEPS = 100       # length of a good race, used to project the evaluation time
_IGNORED_SUFFIXES = (".pyc", ".pyo")


class CheckResult:
    """Outcome of check_submission: errors (the submission will fail), warnings, and measurements."""

    def __init__(self):
        self.errors, self.warnings, self.info = [], [], {}

    @property
    def ok(self):
        return not self.errors

    def __str__(self):
        lines = [f"ERROR: {e}" for e in self.errors] + [f"warning: {w}" for w in self.warnings]
        for k, v in self.info.items():
            lines.append(f"{k}: {v}")
        lines.append("Submission OK." if self.ok else "The submission is NOT valid: fix the errors above.")
        return "\n".join(lines)


def _skip(name):
    return name in _IGNORED or name.startswith(".") or name.endswith(_IGNORED_SUFFIXES)


def make_submission_zip(folder, out_path="submission.zip", check=True):
    """Zip the content of `folder` (agent.py at the root of the zip), without caches and hidden files.
    Returns the CheckResult of the zip if check=True, else None."""
    if not os.path.isfile(os.path.join(folder, AGENT_FILE)):
        raise SubmissionError(f"{folder!r} has no {AGENT_FILE}: give the folder that contains it")
    out_real = os.path.realpath(out_path)
    with zipfile.ZipFile(out_path, "w", zipfile.ZIP_DEFLATED) as z:
        for dirpath, dirnames, filenames in os.walk(folder):
            dirnames[:] = sorted(d for d in dirnames if not _skip(d))
            for name in sorted(filenames):
                path = os.path.join(dirpath, name)
                if _skip(name) or os.path.realpath(path) == out_real:
                    continue
                z.write(path, os.path.relpath(path, folder))
    return check_submission(out_path) if check else None


def _check_structure(submission, result):
    if os.path.isdir(submission):
        names = [os.path.relpath(os.path.join(d, f), submission) for d, _, fs in os.walk(submission) for f in fs]
    elif zipfile.is_zipfile(submission):
        with zipfile.ZipFile(submission) as z:
            names = [n for n in z.namelist() if not n.endswith("/")]
        size_mb = os.path.getsize(submission) / 1e6
        result.info["zip size"] = f"{size_mb:.1f} MB"
        if size_mb > MAX_ZIP_MB:
            result.warnings.append(f"the zip is large ({size_mb:.0f} MB): the upload may fail")
    else:
        result.errors.append(f"{submission!r} is neither a folder nor a zip file")
        return
    if AGENT_FILE not in names:
        nested = [n for n in names if n.endswith("/" + AGENT_FILE) and n.count("/") == 1]
        if nested:
            result.warnings.append(f"{AGENT_FILE} is in the folder {nested[0].split('/')[0]}/ of the zip: it works, "
                                   "but zip the content of the folder (make_submission_zip does it)")
        else:
            result.errors.append(f"no {AGENT_FILE} at the root of the submission")
    junk = sorted({n.split("/")[0] for n in names if any(_skip(p) for p in n.split("/"))})
    if junk:
        result.warnings.append(f"useless files in the submission: {junk[:5]}")


def check_submission(submission, n_steps=TYPICAL_STEPS, seed=None):
    """
    Check a submission (zip or folder) as Codabench will run it: structure, import of agent.py, MyAgent(),
    copy.deepcopy of the agent, valid actions on the first n_steps steps of a training episode, time per step.
    Returns a CheckResult (result.ok, result.errors, result.warnings, result.info); print it to read it.
    """
    result = CheckResult()
    _check_structure(submission, result)
    if result.errors:
        return result
    seed = train_seeds()[0] if seed is None else seed
    with submission_folder(submission) as folder, _working_directory():
        t0 = time.perf_counter()
        try:
            agent, root = load_agent(folder)
        except SubmissionError as exc:
            result.errors.append(str(exc))
            if exc.error is not None and exc.error["type"] in ("FileNotFoundError", "OSError"):
                result.errors.append("the evaluation does not run from your folder: open your files with a path "
                                     "relative to __file__, e.g. os.path.join(os.path.dirname(__file__), 'w.npz')")
            return result
        load_seconds = time.perf_counter() - t0
        result.info["loading time"] = f"{load_seconds:.2f} s"
        if not isinstance(agent, BaseAgent):
            result.warnings.append(f"{AGENT_CLASS} does not derive from sailing_challenge.BaseAgent")
        try:
            player = copy.deepcopy(agent)
        except Exception as exc:
            result.errors.append("every episode plays a fresh copy of the agent, made with copy.deepcopy, which "
                                 f"failed: {format_error(describe_error(exc, root))}. Do not keep open files, "
                                 "generators or locks in the agent's attributes")
            return result
        env = SailingEnv("train")
        obs, info = env.reset(seed=seed)
        times = []
        try:
            player.reset(info["env_info"])
            for step in range(n_steps):
                t = time.perf_counter()
                action = player.act(obs)
                times.append(time.perf_counter() - t)
                a, valid = parse_action(action)
                if not valid:
                    result.errors.append(
                        f"step {step}: act() returned {describe_invalid_action(action)} ({action!r:.60}); expected "
                        f"an int in 0..{N_ACTIONS - 1} ({', '.join(f'{i} = {n}' for i, n in enumerate(ACTION_NAMES))})."
                        " The evaluation plays an invalid action as STAY and counts it")
                    break
                obs, _, terminated, truncated, info = env.step(a)
                if terminated or truncated:
                    break
        except (Exception, SystemExit) as exc:
            result.errors.append(f"the agent raised an exception at step {env.t}: "
                                 f"{format_error(describe_error(exc, root))}")
        if times:
            mean, worst = sum(times) / len(times), max(times)
            result.info["steps checked"] = len(times)
            result.info["time per step"] = f"{1000 * mean:.2f} ms on average, {1000 * worst:.1f} ms at most"
            projected = load_seconds + mean * TYPICAL_STEPS * CODABENCH_EPISODES
            result.info["projected Codabench time"] = (f"about {projected:.0f} s for {CODABENCH_EPISODES} races of "
                                                       f"{TYPICAL_STEPS} steps (budget {BUDGET_SECONDS:.0f} s)")
            if projected > BUDGET_SECONDS:
                result.errors.append(f"too slow: about {projected:.0f} s for the {CODABENCH_EPISODES} episodes of "
                                     f"Codabench, budget {BUDGET_SECONDS:.0f} s. The remaining episodes would score 0")
            elif projected > BUDGET_SECONDS / 2:
                result.warnings.append("more than half of the time budget would be used: little safety margin")
    return result


# ---------------------------------------------------------------------------------------------- command line
def _evaluate_command(args):
    seeds = args.seeds if args.seeds else train_seeds()
    budget = None if args.no_budget else "auto"

    def progress(i, r):
        print(f"  episode {i + 1}/{len(seeds)}: {r['status']}, {r['steps']} steps", file=sys.stderr, flush=True)

    run = evaluate_submission(args.submission, "train", seeds, time_budget=budget,
                              progress=progress if args.verbose else None)
    results, summary = score_run(run, "train", seeds)
    print(f"Train regime, {len(seeds)} seeds. Score: {summary['score']:.2f} (± {summary['score_stderr']:.2f}), "
          f"success rate: {100 * summary['success_rate']:.0f} %, mean arrival time: "
          f"{summary['mean_arrival_time'] if summary['mean_arrival_time'] is None else round(summary['mean_arrival_time'], 2)}"
          f", crashes: {summary['crashes']}.")
    print(format_report(debug_report(results, summary)))
    if run["total_seconds"] and results:
        projected = run["total_seconds"] * CODABENCH_EPISODES / len(results)
        print(f"Projected Codabench time for {CODABENCH_EPISODES} episodes: about {projected:.0f} s "
              f"(budget {BUDGET_SECONDS:.0f} s).")
        if projected > BUDGET_SECONDS:
            print("WARNING: this agent is too slow for the Codabench budget.")
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"summary": summary, "results": results}, f, indent=1)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m sailing_challenge.submission", description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("zip", help="build submission.zip from a folder that contains agent.py, then check it")
    p.add_argument("folder")
    p.add_argument("-o", "--out", default="submission.zip")
    p = sub.add_parser("check", help="check a submission (zip or folder) as Codabench will run it")
    p.add_argument("submission")
    p.add_argument("--steps", type=int, default=TYPICAL_STEPS)
    p = sub.add_parser("evaluate", help="evaluate a submission on the fixed train seeds, with the Codabench code")
    p.add_argument("submission")
    p.add_argument("--seeds", type=int, nargs="*", help="seeds of the training regime (default: train_seeds())")
    p.add_argument("--no-budget", action="store_true", help="no time budget")
    p.add_argument("--json", help="save the results of every episode to this file")
    p.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "zip":
            result = make_submission_zip(args.folder, args.out)
            print(f"Wrote {args.out}")
            print(result)
            return 0 if result.ok else 1
        if args.command == "check":
            result = check_submission(args.submission, n_steps=args.steps)
            print(result)
            return 0 if result.ok else 1
        return _evaluate_command(args)
    except SubmissionError as exc:
        print(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
