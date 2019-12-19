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




v = -70.0
m = m_inf(v)
h = h_inf(v)
n = n_inf(v)
x0 = [v, m, n, h]

if __name__ == "__main__":

    fig, ax = pl.subplots(nrows=3, ncols=2, figsize=(7, 7))
    v = np.arange(-100, 50, 0.01)


    ax[0][0].plot(v, m_inf(v), lw=2, c="k")
    ax[1][0].plot(v, h_inf(v), lw=2, c="k")
    ax[2][0].plot(v, n_inf(v), lw=2, c="k")

    ax[0][1].plot(v, 1.0 / (alpha_m(v) + beta_m(v)), lw=2, c="k")
    # ax[0][1].plot(v, 1.0/(alpha_m1(v)+beta_m(v)), lw=2, c="r")

    ax[1][1].plot(v, 1.0/(alpha_h(v)+beta_h(v)), lw=2, c="k")
    ax[2][1].plot(v, 1.0/(alpha_n(v)+beta_n(v)), lw=2, c="k")

    ax[0][0].set_ylabel(r"$m_{\infty} (v)$")
    ax[1][0].set_ylabel(r"$h_{\infty} (v)$")
    ax[2][0].set_ylabel(r"$n_{\infty} (v)$")

    ax[0][1].set_ylabel(r"$\tau_m [ms]$")
    ax[1][1].set_ylabel(r"$\tau_h [ms]$")
    ax[2][1].set_ylabel(r"$\tau_n [ms]$")

    ax[2][0].set_xlabel("v [mV]", fontsize=14)
    ax[2][1].set_xlabel("v [mV]", fontsize=14)

    for i in range(3):
        for j in range(2):
            ax[i][j].set_xlim(min(v), max(v))
            ax[i][j].set_ylim([0, 1.05])
    ax[1][1].set_ylim([0, 20])
    ax[2][1].set_ylim([0, 20])
    
    
    pl.tight_layout()
    pl.savefig("fig_5_1.png")
    # pl.show()
