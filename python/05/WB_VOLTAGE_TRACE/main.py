""""
Wang-Buzsaki model
Reference:
Gamma Oscillation by Synaptic Inhibition in a Hippocampal 
Interneuronal Network Model
"""

from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl


c = 1.0
g_k = 9.0
g_Na = 35.0
g_l = 0.1
v_k = -90.0
v_na = 55.0
v_l = -65.0
i_ext = 0.75
t_final = 100.0
dt = 0.01
phi = 5.0


def beta_h(v):
    return 1.0 / (exp(-0.1 * (v + 28.0)) + 1.0)

def beta_n(v):
    return 0.125 * exp(-(v + 44.0) / 80.0)

def beta_m(v):
    return 4.0 * exp(-(v + 60.0) / 18.0)

def alpha_h(v):
    return  0.07 * exp(-(v + 58.0) / 20.0)

def alpha_m(v):
    return -0.1 * (v + 35.0) / (exp(-0.1 * (v + 35.0)) - 1.0)

def alpha_n(v):
    return -0.01 * (v + 34.0) / (exp(-0.1 * (v + 34.0)) - 1.0)

def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))

def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))

def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))

def derivative(x0, t):

    v, h, n = x0

    m = alpha_m(v) / (alpha_m(v) + beta_m(v))
    I_Na = g_Na * m ** 3 * h * (v - v_na)
    I_L = g_l * (v - v_l)
    I_K = g_k * n ** 4 * (v - v_k)

    dv = -I_Na - I_K - I_L + i_ext
    dh = phi * (alpha_h(v) * (1-h) - beta_h(v) * h)
    dn = phi * (alpha_n(v) * (1-n) - beta_n(v) * n)

    return [dv, dh, dn]


v = -63.0
# m = m_inf(v)
h = h_inf(v)
n = n_inf(v)
x0 = [v, n, h]


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
    pl.savefig("fig_5_3.png")
    # pl.show()
