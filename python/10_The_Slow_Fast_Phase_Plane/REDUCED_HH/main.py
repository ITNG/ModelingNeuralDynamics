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

    v, n = x0
    
    m = m_inf(v)
    h = 0.83 - n

    I_na = -g_na * h * m ** 3 * (v - v_na)
    I_k = -g_k * n ** 4 * (v - v_k)
    I_l = -g_l * (v - v_l)
    
    dv = (i_ext + I_na + I_k + I_l) / c
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n

    return [dv, dn]


v = -50.0
m = m_inf(v)
n = 0.4
x0 = [v, n]

if __name__ == "__main__":

    fig, ax = pl.subplots(2, figsize=(7, 5), sharex=True)

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]
    n = sol[:, 1]

    ax[0].plot(t, v, lw=2, c="k")
    ax[1].plot(t, n, lw=2, c="r", label="n")
    ax[1].plot(t, 0.83 - n, lw=2, c="g", label="h")
    ax[1].plot(t, m_inf(v), lw=2, c="b", label="m")

    ax[0].set_xlim(min(t), max(t))
    ax[0].set_ylim(-100, 50)
    ax[1].set_ylim(0, 1)
    ax[1].set_xlabel("time [ms]", fontsize=14)
    ax[0].set_ylabel("v [mV]", fontsize=14)
    ax[1].set_ylabel("m, h, n", fontsize=14)
    pl.tight_layout()
    pl.tick_params(labelsize=14)
    pl.savefig("fig_10_2.png")
    # pl.show()
