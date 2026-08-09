from scipy.integrate import odeint
import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 224.
g_na = 112.
g_l = 0.5
v_k = -90.
v_na = 60.
v_l = -70.
i_ext = 7.
t_final = 100.
dt = 0.01


def alpha_h(v):
    return 0.0035 / exp(v / 24.186)


def alpha_m(v):
    return 40 * (75.5 - v) / (exp((75.5 - v) / 13.5) - 1)


def alpha_n(v):
    return (95 - v) / (exp((95 - v) / 11.8) - 1)


def beta_h(v):
    return -0.017 * (v + 51.25) / (exp(-(v + 51.25) / 5.2) - 1)


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


def derivative_3d(x0, t):
    '''three-dimensional model: v, h, n dynamic, m=m_inf(v) quasi-static'''
    v, n, h = x0
    m = m_inf(v)
    dv = (g_k * n ** 2 * (v_k - v) + g_na * m ** 3 * h * (v_na - v)
          + g_l * (v_l - v) + i_ext) / c
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h
    return [dv, dn, dh]


h_plus_n = 0.36


def derivative_2d(x0, t):
    '''reduced (v, n) model: m=m_inf(v), h=h_plus_n-n'''
    v, n = x0
    m = m_inf(v)
    h = h_plus_n - n
    dv = (g_k * n ** 2 * (v_k - v) + g_na * m ** 3 * h * (v_na - v)
          + g_l * (v_l - v) + i_ext) / c
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    return [dv, dn]


v0 = -70.
x0_3d = [v0, n_inf(v0), h_inf(v0)]
x0_2d = [v0, n_inf(v0)]

if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    v_3d = odeint(derivative_3d, x0_3d, t)[:, 0]
    v_2d = odeint(derivative_2d, x0_2d, t)[:, 0]

    fig, ax = plt.subplots(2, figsize=(8, 7), sharex=True)
    ax[0].plot(t, v_3d, lw=2, c='k')
    ax[0].set_ylabel('$v$ [mV]')
    ax[0].set_title('three-dimensional model')

    ax[1].plot(t, v_2d, lw=2, c='k')
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$v$ [mV]')
    ax[1].set_title('two-dimensional model ($h=0.36-n$)')

    for a in ax:
        a.set_xlim(0, t_final)
        a.set_ylim(-95, 55)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
