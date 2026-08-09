import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 36.
g_na = 120.
g_l = 0.3
v_k = -82.
v_na = 45.
v_l = -59.

i_ext = 8.

t_final = 40.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


def alpha_h(v):
    return 0.07 * exp(-(v + 70) / 20)


def alpha_m(v):
    with np.errstate(divide='ignore', invalid='ignore'):
        out = (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))
    return np.where(np.abs(v + 45) > 1e-8, out, 1.0)


def alpha_n(v):
    return 0.01 * (-60. - v) / (exp((-60 - v) / 10) - 1)


def beta_h(v):
    return 1. / (exp(-(v + 40) / 10) + 1)


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def f(v):
    return g_na * m_inf(v) ** 3 * h_inf(v) * (v_na - v) + g_k * n_inf(v) ** 4 * (v_k - v) \
        + g_l * (v_l - v) + i_ext


v_left, v_right = -100., 50.
while v_right - v_left > 1e-12:
    v_c = (v_left + v_right) / 2
    if f(v_c) * f(v_right) < 0:
        v_left = v_c
    else:
        v_right = v_c
v_star = (v_left + v_right) / 2
m_star = m_inf(v_star)
h_star = h_inf(v_star)
n_star = n_inf(v_star)

v = np.empty(m_steps + 1)
m = np.empty(m_steps + 1)
h = np.empty(m_steps + 1)
n = np.empty(m_steps + 1)
v[0], m[0], h[0], n[0] = v_star + 5, m_star, h_star, n_star

for k in range(m_steps):
    v_inc = (g_na * m[k] ** 3 * h[k] * (v_na - v[k]) + g_k * n[k] ** 4 * (v_k - v[k]) + g_l * (v_l - v[k]) + i_ext) / c
    m_inc = alpha_m(v[k]) * (1 - m[k]) - beta_m(v[k]) * m[k]
    h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]
    n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]

    v_tmp = v[k] + dt05 * v_inc
    m_tmp = m[k] + dt05 * m_inc
    h_tmp = h[k] + dt05 * h_inc
    n_tmp = n[k] + dt05 * n_inc

    v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
              + g_l * (v_l - v_tmp) + i_ext) / c
    m_inc = alpha_m(v_tmp) * (1 - m_tmp) - beta_m(v_tmp) * m_tmp
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

    v[k + 1] = v[k] + dt * v_inc
    m[k + 1] = m[k] + dt * m_inc
    h[k + 1] = h[k] + dt * h_inc
    # n is not allowed to drop below n_star -- shows that this decline is
    # what lets the trajectory return to rest
    n[k + 1] = max(n_star, n[k] + dt * n_inc)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    plt.figure(figsize=(7, 3.5))
    plt.plot(t, v, '-k', linewidth=2)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('$v$ [mV]')
    plt.xlim(0, t_final)
    plt.ylim(-100, 50)
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
