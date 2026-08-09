import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 80.
g_na = 100.
g_l = 0.1
v_k = -100.
v_na = 50.
v_l = -67.

g_h = 1.0
v_h = -32.9

t_final = 300.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


def alpha_h(v):
    return 0.128 * exp(-(v + 50) / 18)


def alpha_m(v):
    return 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))


def alpha_n(v):
    return 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))


def beta_h(v):
    return 4. / (1 + exp(-(v + 27) / 5))


def beta_m(v):
    return 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)


def beta_n(v):
    return 0.5 * exp(-(v + 57) / 40)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def r_inf(v):
    return 1. / (1 + exp((v + 84) / 10.2))


def tau_r(v):
    return 1. / (exp(-14.59 - 0.086 * v) + exp(-1.87 + 0.0701 * v))


v = np.zeros(m_steps + 1)
m = np.zeros(m_steps + 1)
h = np.zeros(m_steps + 1)
n = np.zeros(m_steps + 1)
r = np.zeros(m_steps + 1)
v[0], m[0], h[0], n[0], r[0] = -70., 0., 1., 0., 0.5

for k in range(m_steps):
    i_ext = -4. if (t_final / 3 <= k * dt < 2 * t_final / 3) else 0.

    v_inc = (g_k * n[k] ** 4 * (v_k - v[k]) + g_na * m[k] ** 3 * h[k] * (v_na - v[k])
              + g_l * (v_l - v[k]) + g_h * r[k] * (v_h - v[k]) + i_ext) / c
    n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
    h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]
    r_inc = (r_inf(v[k]) - r[k]) / tau_r(v[k])

    v_tmp = v[k] + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h[k] + dt05 * h_inc
    n_tmp = n[k] + dt05 * n_inc
    r_tmp = r[k] + dt05 * r_inc

    v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
              + g_l * (v_l - v_tmp) + g_h * r_tmp * (v_h - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
    r_inc = (r_inf(v_tmp) - r_tmp) / tau_r(v_tmp)

    v[k + 1] = v[k] + dt * v_inc
    m[k + 1] = m_inf(v[k + 1])
    h[k + 1] = h[k] + dt * h_inc
    n[k + 1] = n[k] + dt * n_inc
    r[k + 1] = r[k] + dt * r_inc

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    plt.figure(figsize=(7, 3.5))
    plt.plot(t, v, '-k', linewidth=2)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('$v$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
