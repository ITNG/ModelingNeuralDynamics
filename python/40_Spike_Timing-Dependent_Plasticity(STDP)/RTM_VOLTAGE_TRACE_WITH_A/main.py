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
    v, n, h, a = x0
    dv = i_ext - g_na * h * m_inf(v) ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h
    da = 1 - C * a * (1 + np.tanh(0.1 * v))

    return [dv, dn, dh, da]


v = -70.0
m = m_inf(v)
h = h_inf(v)
n = n_inf(v)
a = 0
x0 = [v, n, h, a]
C = 5

if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]
    a = sol[:, 3]

    fig, ax = pl.subplots(2, figsize=(7, 5), sharex=True)
    ax[0].plot(t, v, lw=2, c="k")
    ax[1].plot(t, a, lw=2, c='k')
    ax[0].set_xlim(min(t), max(t))
    ax[0].set_ylim(-100, 50)
    ax[1].set_xlabel("time [ms]")
    ax[0].set_ylabel("v [mV]")
    ax[1].set_ylabel("a [mV]")
    ax[0].set_yticks(range(-100, 100, 50))
    ax[1].set_ylim(0,20)
    pl.tight_layout()
    pl.savefig("fig_40_3.png")
    # pl.show()
