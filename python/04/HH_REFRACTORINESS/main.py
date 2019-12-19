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


def derivative(x0, t, i_ext, PULSE_ONSET):
    '''
    define HH Model
    '''
    v, m, n, h, = x0

    if (t > PULSE_ONSET) & (t < PULSE_ONSET + 1):
        i_ext = 40.0

    dv = (i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)) / c
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return [dv, dm, dn, dh]



v = -50.0
m = m_inf(v)
h = 0.6 #h_inf(v)
n = 0.4  #n_inf(v)


x0 = [v, m, n, h]

if __name__ == "__main__":

    fig, ax = pl.subplots(3, figsize=(7, 6), sharex=True)

    i_ext = 10.0
    PULSE_ONSETS = [500.0, 5.0, 10.0]
    
    for i in range(3):
        t = np.arange(0, t_final, dt)
        tPulse = PULSE_ONSETS[i]
        sol = odeint(derivative, x0, t, args=(i_ext, tPulse))
        v = sol[:, 0]
        ax[i].plot(t, v, lw=2, c="k")
        ax[i].plot([tPulse, tPulse+1], [-100,-100], lw=3, c="r")
        ax[i].set_ylim([-100, 50])
        ax[i].set_ylabel("v [mV]", fontsize=14)
        ax[i].tick_params(labelsize=14)
    

    ax[0].set_xlim(min(t), max(t))
    ax[2].set_xlabel("time [ms]", fontsize=14)
    
    ax[0].set_yticks(range(-100, 100, 50))

    pl.tight_layout()
    pl.savefig("fig_4_1.png")
    # pl.show()
