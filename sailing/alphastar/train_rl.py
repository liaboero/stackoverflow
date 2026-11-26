"""
Phase 2 of the AlphaStar recipe: reinforcement learning from the imitation policy. The network sails episodes on
fresh random maps (sampling its moves), and is pushed towards the moves that led to better scores than its value
estimate (policy gradient with a learned baseline). As in AlphaStar, a KL penalty keeps it close to the imitation
policy, so that it does not forget what it learned from the expert.

    python alphastar/train_rl.py --init alphastar/bc.pt --out alphastar/rl.pt --iters 40
"""
import argparse, copy, os, sys, time
from multiprocessing import Pool

import numpy as np
import torch
import torch.nn.functional as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge import SailingEnv
from model import PolicyNet, NetAgent
from features import FeatureBuilder

_net = None


def _load(state):
    global _net
    torch.set_num_threads(1)
    if _net is None:
        _net = PolicyNet()
    _net.load_state_dict(state)
    _net.eval()


def rollout(args):
    state, regime, seed = args
    _load(state)
    env = SailingEnv(regime)
    obs, info = env.reset(seed=seed)
    fb = FeatureBuilder(info["env_info"])
    G, V, S, A = [], [], [], []
    done = False
    while not done:
        g, v, s = fb(obs)
        with torch.no_grad():
            logits, _ = _net(torch.from_numpy(g)[None], torch.from_numpy(v)[None], torch.from_numpy(s)[None])
        a = int(torch.distributions.Categorical(logits=logits[0]).sample())
        G.append(g.astype(np.float16)); V.append(v); S.append(s); A.append(a)
        obs, _, term, trunc, _ = env.step(a)
        done = term or trunc
    at = env.arrival_time
    T = len(A)
    # return of each state: the score of the episode seen from that step, 100 * 0.995 ** (arrival - t)
    ret = np.array([100.0 * 0.995 ** (at - t) if at is not None else 0.0 for t in range(T)], dtype=np.float32)
    return np.stack(G), np.stack(V), np.stack(S), np.array(A), ret, (100.0 * 0.995 ** at if at is not None else 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--init", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--iters", type=int, default=40)
    ap.add_argument("--episodes", type=int, default=64)
    ap.add_argument("--kl", type=float, default=0.05)
    ap.add_argument("--lr", type=float, default=1e-4)
    a = ap.parse_args()
    torch.manual_seed(0)
    torch.set_num_threads(os.cpu_count())
    net = PolicyNet()
    net.load_state_dict(torch.load(a.init))
    ref = copy.deepcopy(net).eval()
    for p in ref.parameters():
        p.requires_grad_(False)
    opt = torch.optim.Adam(net.parameters(), lr=a.lr)
    seed = 300000   # disjoint from the demonstrations (200000+) and from the test seeds
    with Pool(os.cpu_count()) as pool:
        for it in range(a.iters):
            t0 = time.time()
            state = {k: v.clone() for k, v in net.state_dict().items()}
            jobs = [(state, "train" if (seed + i) % 2 == 0 else "random", seed + i) for i in range(a.episodes)]
            seed += a.episodes
            res = pool.map(rollout, jobs, chunksize=1)
            G = torch.from_numpy(np.concatenate([r[0] for r in res])).float()
            V = torch.from_numpy(np.concatenate([r[1] for r in res]))
            S = torch.from_numpy(np.concatenate([r[2] for r in res]))
            A = torch.from_numpy(np.concatenate([r[3] for r in res]))
            R = torch.from_numpy(np.concatenate([r[4] for r in res])) / 100.0
            scores = [r[5] for r in res]
            net.train()
            for epoch in range(2):
                perm = torch.randperm(len(A))
                for i in range(0, len(A), 512):
                    idx = perm[i:i + 512]
                    logits, value = net(G[idx], V[idx], S[idx])
                    adv = (R[idx] - value).detach()
                    adv = (adv - adv.mean()) / (adv.std() + 1e-6)
                    logp = F.log_softmax(logits, 1)
                    pg = -(logp.gather(1, A[idx][:, None]).squeeze(1) * adv).mean()
                    vloss = F.mse_loss(value, R[idx])
                    with torch.no_grad():
                        ref_logits, _ = ref(G[idx], V[idx], S[idx])
                    kl = (logp.exp() * (logp - F.log_softmax(ref_logits, 1))).sum(1).mean()
                    loss = pg + vloss + a.kl * kl
                    opt.zero_grad(); loss.backward()
                    torch.nn.utils.clip_grad_norm_(net.parameters(), 1.0)
                    opt.step()
            torch.save(net.state_dict(), a.out)
            print(f"iter {it+1}  sampled-policy score {np.mean(scores):.2f}  kl {float(kl):.3f}  vloss {float(vloss):.4f}  {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
