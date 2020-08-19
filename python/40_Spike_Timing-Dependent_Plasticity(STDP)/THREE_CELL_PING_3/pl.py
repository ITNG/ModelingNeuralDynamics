import numpy as np
import pylab as pl
from lib import read_from_file
from main import num_e, num_i, t_final


data = np.loadtxt("data.txt")
g = data[:, 0]
delta = data[:, 1]
f = data[:, 2]


fig, ax = pl.subplots(ncols=2, figsize=(8, 3.5))

ax[0].plot(g, delta, lw=2, c='k')
ax[1].plot(g, f, lw=2, c='k')

ax[0].set_ylabel(r"$\Delta$", fontsize=13)
ax[1].set_ylabel(r"$f$", fontsize=13)
for i in range(2):

    ax[i].set_xlabel(r"$\bar{g_{EE}}$", fontsize=13)
    ax[i].tick_params(labelsize=13)

pl.margins(x=0.01)
pl.tight_layout()
pl.savefig("fig.png")


    
