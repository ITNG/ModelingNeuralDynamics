from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from model import simulate_poisson_ping


def simulate(seed=63806, t_final=500.0):
    return simulate_poisson_ping(seed=seed, t_final=t_final)


def plot(result, t_final=500.0):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(result["t_i_spikes"], result["i_i_spikes"], ".b", markersize=2)
    ax.plot(result["t_e_spikes"], result["i_e_spikes"] + 50, ".r", markersize=2)
    ax.axhline(50.5, color="k", linestyle="--", linewidth=1)
    t = np.linspace(0, t_final, 1001)
    ax.plot(t, np.exp(4 * np.cos(np.pi * t / 31) ** 2) - 1 + 60, "-g", linewidth=2)
    ax.set(xlim=(0, t_final), ylim=(0, 251), xlabel="$t$ [ms]", yticks=[50, 250])
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    result = simulate()
    print(f"f_hat_e = {result['f_hat_e']:.3f}")
    print(f"f_hat_i = {result['f_hat_i']:.3f}")
    plot(result).savefig("fig.png", dpi=150)
