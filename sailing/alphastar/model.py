"""Policy network (AlphaStar-style encoder, scaled down): a CNN over the wind maps + an MLP over the boat state."""
import os, sys

import numpy as np
import torch
import torch.nn as nn

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.dirname(HERE), HERE]
from sailing_challenge.agents import BaseAgent
from features import FeatureBuilder

C_IN, V_IN = 11, 44


class PolicyNet(nn.Module):
    def __init__(self, width=48):
        super().__init__()
        self.cnn = nn.Sequential(
            nn.Conv2d(C_IN, width, 3, padding=1), nn.ReLU(),
            nn.Conv2d(width, width, 3, stride=2, padding=1), nn.ReLU(),      # 16
            nn.Conv2d(width, 2 * width, 3, stride=2, padding=1), nn.ReLU(),  # 8
            nn.Conv2d(2 * width, 2 * width, 3, stride=2, padding=1), nn.ReLU(),  # 4
            nn.Flatten(), nn.Linear(2 * width * 16, 256), nn.ReLU())
        self.mlp = nn.Sequential(nn.Linear(V_IN, 128), nn.ReLU(), nn.Linear(128, 128), nn.ReLU())
        self.head = nn.Sequential(nn.Linear(384, 256), nn.ReLU(), nn.Linear(256, 9))
        self.value = nn.Sequential(nn.Linear(384, 128), nn.ReLU(), nn.Linear(128, 1))

    def forward(self, grid, vec, safe=None):
        z = torch.cat([self.cnn(grid), self.mlp(vec)], 1)
        logits = self.head(z)
        if safe is not None:   # never choose a move that hits an island
            logits = logits.masked_fill(~safe, -1e9)
        return logits, self.value(z).squeeze(1)


class NetAgent(BaseAgent):
    """The trained network as an agent: greedy action (or sampled, for reinforcement learning)."""

    def __init__(self, path_or_net, sample=False):
        torch.set_num_threads(1)
        if isinstance(path_or_net, str):
            self.net = PolicyNet()
            self.net.load_state_dict(torch.load(path_or_net, map_location="cpu"))
        else:
            self.net = path_or_net
        self.net.eval()
        self.sample = sample

    def reset(self, env_info):
        self.fb = FeatureBuilder(env_info)

    def act(self, obs):
        g, v, s = self.fb(obs)
        with torch.no_grad():
            logits, _ = self.net(torch.from_numpy(g)[None], torch.from_numpy(v)[None], torch.from_numpy(s)[None])
        if self.sample:
            return int(torch.distributions.Categorical(logits=logits[0]).sample())
        return int(logits[0].argmax())
