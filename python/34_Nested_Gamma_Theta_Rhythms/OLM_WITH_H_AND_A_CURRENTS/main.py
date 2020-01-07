from scipy.integrate import odeint
import numpy as np
import pylab as pl
from numpy import exp
import pylab as pl
from lib import *


c = 1.3
g_k = 23.0
g_na = 30.0
g_h = 12.0
g_A = 22.0
g_l = 0.05
v_k = -100.0
v_na = 90.0
v_l = -70.0
v_h = -32.9
v_A = -90.0

i_ext = 0

t_final = 500.0
dt = 0.01


if __name__ == "__main__":

    v = -63.0
    x0 = [v, h_o_inf(v), n_o_inf(v), r_o_inf(v), a_o_inf(v), b_o_inf(v)]
    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]
    a = sol[:, 4]
    b = sol[:, 5]

    fig, ax = pl.subplots(2, figsize=(7, 4), sharex=True)

    ax[0].plot(t, v, lw=2, c="k")
    ax[0].set_xlim(min(t), max(t))
    ax[0].set_ylabel("v [mV]")
    ax[0].set_yticks(range(-100, 100, 50))

    ax[1].plot(t, a * b, lw=2, c="k")
    ax[1].set_ylim(0, 0.05)
    ax[1].set_xlabel("time [ms]", fontsize=14)
    ax[1].set_ylabel("ab", fontsize=14)

    for i in range(2):
        ax[i].tick_params(labelsize=14)

    pl.tight_layout()
    pl.savefig("fig_34_8.png")
    pl.close()
