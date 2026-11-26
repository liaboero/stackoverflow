"""
Write replay files for the replay page, from submission folders (agent.py with a class MyAgent):

    python -m sailing_challenge.make_replays my_submission/ [other_submission/ ...] --seeds 1-3 --out replays

One file per seed, out/replay_seed<seed>.json, where all the agents race. Then open replay/index.html in a browser
and choose the files. From a notebook, use sailing_challenge.export_replays(...) with agent objects instead.
"""
import argparse
import os

from .evaluation import load_agent
from .replay import export_replays


def parse_seeds(text):
    """'1-3,7' -> [1, 2, 3, 7]."""
    seeds = []
    for part in text.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-", 1)
            seeds.extend(range(int(a), int(b) + 1))
        elif part:
            seeds.append(int(part))
    return seeds


def load_submission_agent(folder):
    """MyAgent() from folder/agent.py, loaded exactly like the evaluation does."""
    return load_agent(folder)[0]


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m sailing_challenge.make_replays", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("submissions", nargs="+", help="submission folders (each with agent.py and a class MyAgent)")
    parser.add_argument("--seeds", default="1-5", help="seeds of the training regime, e.g. 1-5 or 1,4,9 (default 1-5)")
    parser.add_argument("--out", default="replays", help="output folder (default: replays)")
    parser.add_argument("--names", help="names of the boats, comma separated (default: the folder names)")
    parser.add_argument("--static-wind", action="store_true", help="only the base wind, which never changes")
    args = parser.parse_args(argv)
    names = args.names.split(",") if args.names else [os.path.basename(os.path.normpath(s)) for s in args.submissions]
    if len(names) != len(args.submissions):
        parser.error("one name per submission")
    if len(set(names)) != len(names):
        names = [f"{n} ({k + 1})" for k, n in enumerate(names)]
    agents = {name: load_submission_agent(folder) for name, folder in zip(names, args.submissions)}
    export_replays(agents, seeds=parse_seeds(args.seeds), out_dir=args.out, static_wind=args.static_wind)


if __name__ == "__main__":
    main()
