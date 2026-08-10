from pathlib import Path
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from model import simulate_two_cell


def simulate(t_final=200.0):
    coupled = simulate_two_cell(g_ei=0.5, t_final=t_final)
    mean_s_i = float(coupled["s_i"][1:].mean())
    averaged = simulate_two_cell(g_ei=0.5, t_final=t_final, mean_inhibition=mean_s_i)
    return {"coupled": coupled, "mean_inhibition": averaged, "mean_s_i": mean_s_i}


def plot(result):
    run = result["mean_inhibition"]
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(run["t"], run["i_main"], "-k", linewidth=2)
    axes[0].plot(run["t"], run["i_dist"], "--k", linewidth=2)
    axes[0].set_ylim(0, 4)
    axes[1].plot(run["t_e_spikes"], 2 + 0 * run["t_e_spikes"], ".r", markersize=14)
    axes[1].plot(run["t_i_spikes"], 1 + 0 * run["t_i_spikes"], ".b", markersize=14)
    for time in range(0, round(run["t"][-1]) + 1, 25):
        axes[1].axvline(time, color="k", linestyle="--", linewidth=1)
    axes[1].set(xlim=(0, run["t"][-1]), ylim=(0, 3), xlabel="$t$ [ms]", yticks=[])
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    plot(simulate()).savefig("fig.png", dpi=150)
