import numpy as np
import pylab as pl
from lib import read_from_file
from main import num_e, num_i

t_e_spikes = read_from_file("t_e_spikes.txt")
t_i_spikes = read_from_file("t_i_spikes.txt")

fig, ax = pl.subplots(1, figsize=(7, 2))

for i in range(num_i):
    ax.plot(t_i_spikes[i], [i] * len(t_i_spikes[i]), "b.")

for i in range(num_i, num_i + num_e):
    ax.plot(t_e_spikes[i - num_i], [i] *
            len(t_e_spikes[i - num_i]), "r.")


ax.set_xlabel("time [ms]")
ax.set_ylabel("neuron #")

pl.tight_layout()
pl.savefig("fig.png")

