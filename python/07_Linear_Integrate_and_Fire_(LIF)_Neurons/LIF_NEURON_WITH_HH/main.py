from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl

c = 1.0
g_k = 36.0
g_na = 120.0
g_l = 0.3
v_k = -82.0
v_na = 45.0
v_l = -59.0
i_ext = 7.0
t_final = 100.0
dt = 0.01


def alpha_h(v):
    return 0.07 * exp(-(v + 70) / 20)


def alpha_m(v):
    return (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))


def alpha_n(v):
    return 0.01 * (-60.0 - v) / (exp((-60 - v) / 10) - 1)


def beta_h(v):
    return 1. / (exp(-(v + 40) / 10) + 1)


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def derivative(x0, t):

    v, m, n, h, = x0
    I_na = -g_na * h * m ** 3 * (v - v_na)
    I_k = -g_k * n ** 4 * (v - v_k)
    I_l = -g_l * (v - v_l)
    
    dv = (i_ext + I_na + I_k + I_l) / c
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return [dv, dm, dn, dh]


v = -20.0
m = m_inf(v)
h = h_inf(v)
n = n_inf(v)
x0 = [v, m, n, h]

if __name__ == "__main__":

    pl.figure(figsize=(7, 3))

    v = np.asarray([-82.0, -54.0])
    t = np.asarray([0.0, (v[1] - v[0]) / 1.65])
    pl.plot(t, v, c="b", lw=2)

    for i in range(1, 6):
        pl.plot(t + t[1] * i, v, c="b", lw=2)
        pl.plot([(i)*t[1], (i)*t[1]], v, c="b", lw=2, ls="--")



    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]

    pl.plot(t, v, lw=2, c="k")
    pl.xlim(min(t), max(t))
    pl.ylim(-100, 50)
    pl.xlabel("time [ms]")
    pl.ylabel("v [mV]")
    pl.yticks(range(-100, 100, 50))
    pl.tight_layout()
    pl.savefig("fig_7_1.png")
    # pl.show()
