"""Interface of an agent. The evaluation calls reset(env_info) at the start of every episode, then act(obs)."""


class BaseAgent:
    def reset(self, env_info):
        """Start of an episode. env_info holds the fixed characteristics of the environment (islands, physics...)."""

    def act(self, obs):
        """Return an action in 0..8 (0 = N, 1 = NE, ..., 7 = NW, 8 = stay)."""
        raise NotImplementedError

    def seed(self, seed=None):
        """Optional: seed the agent's own randomness (never the seed of the episode)."""
