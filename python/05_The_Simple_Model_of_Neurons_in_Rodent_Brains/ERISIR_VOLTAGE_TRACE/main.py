from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl

c = 1
g_k = 224.0
g_na = 112
g_l = 0.5
v_k = -90.0
v_na = 60
v_l = -70
i_ext = 7.0
t_final = 100
dt = 0.01


def alpha_h(v):
    return 0.0035 / exp(v / 24.186)
    


def alpha_m(v):
    return 40 * (75.5 - v) / (exp((75.5 - v) / 13.5) - 1)
    


def alpha_n(v):
    return (95 - v) / (exp((95 - v) / 11.8) - 1)


def beta_h(v):
    return - 0.017 * (v + 51.25) / (exp(-(v + 51.25) / 5.2) - 1)


def beta_m(v):
    return 1.2262 / exp(v / 42.248)
    


def beta_n(v):
    return 0.025 / exp(v / 22.222)


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
    v, n, h, = x0
    m = alpha_m(v) / (alpha_m(v) + beta_m(v))
    dv = i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 2 * (v - v_k) - g_l * (v - v_l)
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return [dv, dn, dh]


v = -70.0
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
    pl.savefig("fig_5_4.png")
    # pl.show()
