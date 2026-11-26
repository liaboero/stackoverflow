"""
Phase 1 of the AlphaStar recipe: supervised imitation of expert games. Here the expert is the Oracle (it knows the
future wind), the student only sees the present and the past: to imitate it, it has to learn to anticipate the wind.

    python alphastar/train_bc.py --demos demos.npz --out alphastar/bc.pt
"""
import argparse, os, sys, time

import numpy as np
import torch
import torch.nn.functional as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from model import PolicyNet


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--demos", required=True, nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--init", default=None)
    a = ap.parse_args()
    torch.manual_seed(0)
    torch.set_num_threads(os.cpu_count())
    parts = [np.load(p) for p in a.demos]
    grid = np.concatenate([d["grid"] for d in parts])
    vec = np.concatenate([d["vec"] for d in parts])
    safe = np.concatenate([d["safe"] for d in parts])
    act = np.concatenate([d["act"] for d in parts])
    ep = np.concatenate([d["ep"] + 100000 * i for i, d in enumerate(parts)])
    eps = np.unique(ep)
    rng = np.random.default_rng(0)
    val_eps = set(rng.choice(eps, max(1, len(eps) // 10), replace=False).tolist())
    is_val = np.array([e in val_eps for e in ep])
    tr, va = np.where(~is_val)[0], np.where(is_val)[0]
    print(f"samples {len(act)}  train {len(tr)}  val {len(va)}", flush=True)
    G, V, S, A = (torch.from_numpy(x) for x in (grid, vec, safe, act))
    net = PolicyNet()
    if a.init:
        net.load_state_dict(torch.load(a.init))
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, a.epochs)
    bs = 256
    for epoch in range(a.epochs):
        net.train()
        t0, perm, tot, n = time.time(), rng.permutation(tr), 0.0, 0
        for i in range(0, len(perm), bs):
            idx = torch.from_numpy(perm[i:i + bs])
            logits, _ = net(G[idx].float(), V[idx], S[idx])
            loss = F.cross_entropy(logits, A[idx])   # masked logits: no label smoothing (it would push mass onto forbidden moves)
            opt.zero_grad(); loss.backward(); opt.step()
            tot += float(loss.detach()) * len(idx); n += len(idx)
        sched.step()
        net.eval()
        with torch.no_grad():
            correct = 0
            for i in range(0, len(va), 1024):
                idx = torch.from_numpy(va[i:i + 1024])
                logits, _ = net(G[idx].float(), V[idx], S[idx])
                correct += int((logits.argmax(1) == A[idx]).sum())
        print(f"epoch {epoch+1}  loss {tot/n:.3f}  val-acc {correct/len(va):.3f}  {time.time()-t0:.0f}s", flush=True)
        torch.save(net.state_dict(), a.out)


if __name__ == "__main__":
    main()
