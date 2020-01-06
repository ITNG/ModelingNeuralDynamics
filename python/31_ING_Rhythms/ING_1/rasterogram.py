import numpy as np
import pylab as pl
from lib import read_from_file
from main import num_i

t_i_spikes = read_from_file("t_i_spikes.txt")
lfp = np.loadtxt("lfp.txt", dtype=float)

fig, ax = pl.subplots(2, figsize=(10, 8), sharex=True)

for i in range(num_i):
    ax[0].plot(t_i_spikes[i], [i] * len(t_i_spikes[i]), "b.")

t = lfp[:, 0]
lfp = lfp[:, 1]

ax[1].plot(t, lfp, lw=2, color="k")

ax[1].set_xlabel("time [ms]", fontsize=18)
ax[0].set_ylabel("neuron #", fontsize=18)
ax[1].set_ylabel("mean(v), I-cells", fontsize=18)

for i in range(2):
    ax[i].tick_params(labelsize=14)
ax[1].set_ylim(-100, 50)
ax[1].set_xlim(0, np.max(t))
pl.tight_layout()
pl.savefig("fig.png")
pl.show()
