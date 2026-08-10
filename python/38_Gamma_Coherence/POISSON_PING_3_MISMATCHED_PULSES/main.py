from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from model import simulate_poisson_ping

P = 29.0


def simulate_phases(seed=63806, phases=np.arange(0.1, 1.0, 0.1), t_final=500.0):
    counts = []
    for index, phase in enumerate(phases):
        result = simulate_poisson_ping(seed=seed + index, t_final=t_final, pulse_period=P, pulse_phase=phase)
        counts.append(np.count_nonzero(result["i_e_spikes"] < 5))
    counts = np.asarray(counts)
    slope, intercept = np.polyfit(phases, counts, 1)
    return {"phases": np.asarray(phases), "spike_counts": counts, "fit_slope": slope, "fit_intercept": intercept}


def plot(result):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(result["phases"], result["spike_counts"], "ok", markerfacecolor="none", markersize=8)
    ax.plot(result["phases"], result["fit_slope"] * result["phases"] + result["fit_intercept"], "-g", linewidth=3)
    ax.set(xlim=(0, 1), ylim=(10, max(40, result["spike_counts"].max() + 2)), xlabel=r"$\varphi$", ylabel="# spikes in E-cells 1-5")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    plot(simulate_phases()).savefig("fig.png", dpi=150)
