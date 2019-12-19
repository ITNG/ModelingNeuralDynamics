from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl

c = 1
g_k = 80
g_na = 100
g_l = 0.1
v_k = -100
v_na = 50
v_l = -67
i_ext = 1.5
t_final = 100
dt = 0.01


def alpha_h(v):
    return 0.128 * exp(-(v + 50.0) / 18.0)


def alpha_m(v):
    return 0.32 * (v + 54) / (1.0 - exp(-(v + 54.0) / 4.0))


def alpha_n(v):
    return 0.032 * (v + 52) / (1.0 - exp(-(v + 52.0) / 5.0))


def beta_h(v):
    return 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))


def beta_m(v):
    return 0.28 * (v + 27.0) / (exp((v + 27.0) / 5.0) - 1.0)


def beta_n(v):
    return 0.5 * exp(-(v + 57.0) / 40.0)


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def derivative(x0, t):
    '''
    define Traub Model
    '''
    v, m, n, h, = x0
    dv = i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return [dv, dm, dn, dh]


v = -70.0
m = m_inf(v)
h = h_inf(v)
n = n_inf(v)
x0 = [v, m, n, h]

if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]

    pl.figure(figsize=(7, 3))
    pl.plot(t, v, lw=2, c="k")
    pl.xlim(min(t), max(t))
    pl.ylim(-100, 50)
    pl.xlabel("time [ms]")
    pl.ylabel("v [mV]")
    pl.yticks(range(-100, 100, 50))
    pl.tight_layout()
    pl.savefig("fig_5_2.png")
    # pl.show()
