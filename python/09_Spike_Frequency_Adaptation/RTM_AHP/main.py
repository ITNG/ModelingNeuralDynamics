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
g_ahp = 0.25
i_ext = 1.5
t_final = 300
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


def cac_inf(v):
    return (120 - v) / (1 + exp(-(v + 15) / 5)) * 4 / 25


def derivative(x0, t):
    '''
    define Traub Model with a calcium-dependent AHP current

    m is not its own state -- matlab/09/RTM_AHP/make_figure.m always sets
    m = m_inf(v) (quasi-static), never integrates it as an ODE. The [Ca2+]
    time constant is a plain constant (80), not the voltage-dependent
    tau_w -- matlab/09/RTM_AHP/tau_w.m exists but is unused there.
    '''
    v, n, h, ca = x0
    m = m_inf(v)
    dv = i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l) \
        - g_ahp * ca * (v - v_k)
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h
    dca = (cac_inf(v) - ca) / 80.0

    return [dv, dn, dh, dca]


v = -70.0
h = h_inf(v)
n = n_inf(v)
ca = 0.0
x0 = [v, n, h, ca]

if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]
    ca = sol[:, 3]

    fig, ax = pl.subplots(2, figsize=(7, 6), sharex=True)
    ax[0].plot(t, v, lw=2, c="k")
    ax[0].set_ylabel("v [mV]")

    ax[1].plot(t, ca, lw=2, c="k")
    ax[1].set_xlabel("t [ms]")
    ax[1].set_ylabel("$[Ca^{2+}]$")
    ax[1].set_ylim(0, max(ca) * 1.2)

    ax[0].set_xlim(min(t), max(t))
    pl.tight_layout()
    pl.savefig("fig.png")
    # pl.show()
