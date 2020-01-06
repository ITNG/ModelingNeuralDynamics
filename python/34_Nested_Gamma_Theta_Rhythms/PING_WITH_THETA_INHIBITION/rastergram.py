import numpy as np
import pylab as pl
from lib import read_from_file
from main import num_e, num_i

t_e_spikes = read_from_file("t_e_spikes.txt")
t_i_spikes = read_from_file("t_i_spikes.txt")
lfp_v = np.loadtxt("lfp_v.txt", dtype=float)
lfp_s = np.loadtxt("lfp_s.txt", dtype=float)

fig, ax = pl.subplots(3, figsize=(10, 8), sharex=True)

for i in range(num_i):
    ax[0].plot(t_i_spikes[i], [i] * len(t_i_spikes[i]), "b.")

for i in range(num_i, num_i + num_e):
    ax[0].plot(t_e_spikes[i - num_i], [i] *
            len(t_e_spikes[i - num_i]), "r.")
t = lfp_v[:, 0]
lfp_v = lfp_v[:, 1]

ax[1].plot(t, lfp_v, lw=2, color="k")

ax[2].set_xlabel("time [ms]", fontsize=18)
ax[0].set_ylabel("neuron #", fontsize=18)
ax[1].set_ylabel(r"mean($v_E$)", fontsize=18)


t = lfp_s[:, 0]
lfp_s = lfp_s[:, 1]

ax[2].plot(t, lfp_s, lw=2, color="k")
ax[2].set_ylabel(r"mean($s_E$)", fontsize=18)

for i in range(3):
    ax[i].tick_params(labelsize=14)
# ax[1].set_ylim(-100, 50)
ax[0].set_xlim(0, np.max(t))
pl.tight_layout()
pl.savefig("fig.png")
# pl.show()

