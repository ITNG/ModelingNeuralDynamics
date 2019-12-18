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
i_ext = 10.0
t_final = 50.0
dt = 0.01


def beta_n(v):
    return 0.125 * exp(-(v + 70.0) / 80.0)


def beta_m(v):
    return 4.0 * exp(-(v + 70.0) / 18.0)


def beta_h(v):
    return 1. / (exp(-(v + 40.0) / 10.0) + 1.0)


def alpha_n(v):
    return 0.01 * (-60.0 - v) / (exp((-60.0 - v) / 10.0) - 1.0)


def alpha_m(v):
    if np.abs(v+45.0) > 1.0e-8:
        return (v + 45.0) / 10.0 / (1.0 - exp(-(v + 45.0) / 10.0))
    else:
        return 1.0


def alpha_h(v):
    return 0.07*exp(-(v+70)/20)


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def derivative(x0, t):
    '''
    define HH Model
    '''
    v, m, n, h, = x0
    dv = (i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)) / c
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return [dv, dm, dn, dh]



v = -50.0
m = m_inf(v)
h = 0.6 #h_inf(v)
n = 0.4 #n_inf(v)

x0 = [v, m, n, h]

if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]
    m = sol[:, 1]
    n = sol[:, 2]
    h = sol[:, 3]

    fig, ax = pl.subplots(1, figsize=(5, 5))
    
    ax.plot(v, n, lw=2, c="k")

    ax.set_xlim(-100, 50)
    ax.set_ylim([0, 1])

    ax.set_xlabel("v [mV]")
    ax.set_ylabel("n")

    pl.savefig("fig_4_2.png")
