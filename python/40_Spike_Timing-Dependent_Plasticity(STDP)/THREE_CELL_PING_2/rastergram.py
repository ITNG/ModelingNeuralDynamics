import numpy as np
import pylab as pl
from lib import read_from_file
from main import num_e, num_i


def plot_raster(e_file, i_file, ax):
    t_e_spikes = read_from_file(e_file)
    t_i_spikes = read_from_file(i_file)

    for i in range(num_i):
        if len(t_i_spikes[i]) > 0:
            ax.plot(t_i_spikes[i], [i] * len(t_i_spikes[i]), "b.")

    for i in range(num_i, num_i + num_e):
        if len(t_e_spikes) > 0:
            ax.plot(t_e_spikes[i - num_i], [i] *
                    len(t_e_spikes[i - num_i]), "r.")


    ax.set_xlabel("time [ms]")
    ax.set_ylabel("neuron #")

fig, ax = pl.subplots(2, figsize=(7, 3), sharex=True)
plot_raster('t_e_spikes1.txt', 't_i_spikes1.txt', ax[0])
plot_raster('t_e_spikes2.txt', 't_i_spikes2.txt', ax[1])
pl.margins(x=0.01)
pl.tight_layout()
pl.savefig("fig.png")
