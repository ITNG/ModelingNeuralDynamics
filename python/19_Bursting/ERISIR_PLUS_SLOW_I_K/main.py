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

i_ext = 7.5

g_k_slow = 1.5
tau_n_slow = 100.

t_final = 1000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


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


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def n_slow_inf(v):
    return 1. / (1 + exp((-20 - v) / 5))


v = np.zeros(m_steps + 1)
m = np.zeros(m_steps + 1)
h = np.zeros(m_steps + 1)
n = np.zeros(m_steps + 1)
n_slow = np.zeros(m_steps + 1)
v[0] = -70.
m[0] = m_inf(v[0])
h[0] = h_inf(v[0])
n[0] = n_inf(v[0])
n_slow[0] = n_slow_inf(v[0])

for k in range(m_steps):
    v_inc = (g_k * n[k] ** 2 * (v_k - v[k]) + g_na * m[k] ** 3 * h[k] * (v_na - v[k])
              + g_k_slow * n_slow[k] * (v_k - v[k]) + g_l * (v_l - v[k]) + i_ext) / c
    n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
    h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]
    n_slow_inc = (n_slow_inf(v[k]) - n_slow[k]) / tau_n_slow

    v_tmp = v[k] + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h[k] + dt05 * h_inc
    n_tmp = n[k] + dt05 * n_inc
    n_slow_tmp = n_slow[k] + dt05 * n_slow_inc

    v_inc = (g_k * n_tmp ** 2 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
              + g_k_slow * n_slow_tmp * (v_k - v_tmp) + g_l * (v_l - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
    n_slow_inc = (n_slow_inf(v_tmp) - n_slow_tmp) / tau_n_slow

    v[k + 1] = v[k] + dt * v_inc
    m[k + 1] = m_inf(v[k + 1])
    h[k + 1] = h[k] + dt * h_inc
    n[k + 1] = n[k] + dt * n_inc
    n_slow[k + 1] = n_slow[k] + dt * n_slow_inc

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    plt.figure(figsize=(7, 4))
    plt.plot(t, v, '-k', linewidth=2)
    plt.xlim(0, t_final)
    plt.ylim(-95, 55)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('$v$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
