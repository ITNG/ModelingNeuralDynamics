from scipy.integrate import odeint
import numpy as np
import pylab as pl
from numpy import exp
import pylab as pl
from lib import *

c = 1.3
g_k = 23
g_na = 30
g_l = 0.05
v_k = -100
v_na = 90
v_l = -70

i_ext = 1.5
t_final = 200.0
dt = 0.01



if __name__ == "__main__":

    v = -63.0
    x0 = [v, h_inf(v), n_inf(v)]
    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]

    pl.figure(figsize=(7, 3))
    pl.plot(t, v, lw=2, c="k")
    pl.xlim(min(t), max(t))
    # pl.ylim(-100, 50)
    pl.xlabel("time [ms]")
    pl.ylabel("v [mV]")
    pl.yticks(range(-100, 100, 50))
    pl.tight_layout()
    pl.savefig("fig_34_5.png")
    pl.close()
    