from scipy.integrate import odeint
import numpy as np
import pylab as pl
from numpy import exp
import pylab as pl
from lib import *

c = 1.3
g_k = 23
g_na = 30
g_h = 12
g_l = 0.05
v_k = -100
v_na = 90
v_l = -70
v_h = -32.9

i_ext = 0.0
t_final = 200.0
dt = 0.01


if __name__ == "__main__":

    v = -63.0
    x0 = [v, h_inf(v), n_inf(v), r_inf(v)]
    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]
    r = sol[:, -1]

    fig, ax = pl.subplots(2, figsize=(7, 4), sharex=True)

    ax[0].plot(t, v, lw=2, c="k")
    ax[0].set_xlim(min(t), max(t))
    # pl.ylim(-100, 50)
    ax[0].set_ylabel("v [mV]")
    ax[0].set_yticks(range(-100, 100, 50))

    ax[1].plot(t, r, lw=2, c="k")
    ax[1].set_ylim(0, 0.01)
    ax[1].set_xlabel("time [ms]")
    ax[1].set_ylabel("r [mV]")
    ax[1].set_yticks([0., 0.005, 0.01])
    

    pl.tight_layout()
    pl.savefig("fig_34_5.png")
    pl.close()
