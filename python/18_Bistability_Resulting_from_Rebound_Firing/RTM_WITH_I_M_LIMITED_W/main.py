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
g_m = 0.2

i_ext = 0.506

t_final = 500.
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


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def w_inf(v):
    return 1. / (1 + exp(-(v + 35) / 10))


def tau_w(v):
    return 400. / (3.3 * exp((v + 35) / 20) + exp(-(v + 35) / 20))


def f(v):
    return g_na * m_inf(v) ** 3 * h_inf(v) * (v_na - v) + g_k * n_inf(v) ** 4 * (v_k - v) \
        + g_l * (v_l - v) + g_m * w_inf(v) * (v_k - v) + i_ext


v_grid = -100 + np.arange(100001) / 100000 * 150
with np.errstate(invalid='ignore'):
    f_grid = f(v_grid)
ind = np.where(f_grid[:-1] * f_grid[1:] <= 0)[0].min()
v_left, v_right = v_grid[ind], v_grid[ind + 1]
while v_right - v_left > 1e-14:
    v_c = (v_left + v_right) / 2
    if f(v_c) * f(v_right) <= 0:
        v_left = v_c
    else:
        v_right = v_c
v_star = (v_left + v_right) / 2
m_star = m_inf(v_star)
h_star = h_inf(v_star)
n_star = n_inf(v_star)
w_star = w_inf(v_star)

v = np.empty(m_steps + 1)
m = np.empty(m_steps + 1)
h = np.empty(m_steps + 1)
n = np.empty(m_steps + 1)
w = np.empty(m_steps + 1)
v[0], m[0], h[0], n[0], w[0] = v_star + 1, m_star, h_star, n_star, w_star

for k in range(m_steps):
    v_inc = (g_na * m[k] ** 3 * h[k] * (v_na - v[k]) + g_k * n[k] ** 4 * (v_k - v[k]) + g_l * (v_l - v[k])
              + g_m * w[k] * (v_k - v[k]) + i_ext) / c
    h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]
    n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
    w_inc = (w_inf(v[k]) - w[k]) / tau_w(v[k])

    v_tmp = v[k] + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h[k] + dt05 * h_inc
    n_tmp = n[k] + dt05 * n_inc
    w_tmp = w[k] + dt05 * w_inc

    v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
              + g_l * (v_l - v_tmp) + g_m * w_tmp * (v_k - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
    w_inc = (w_inf(v_tmp) - w[k]) / tau_w(v_tmp)

    v[k + 1] = v[k] + dt * v_inc
    m[k + 1] = m_inf(v[k + 1])
    h[k + 1] = h[k] + dt * h_inc
    n[k + 1] = n[k] + dt * n_inc
    # w is not allowed to drop below w_star -- shows that this decline is
    # what lets the trajectory return to repetitive firing
    w[k + 1] = max(w[k] + dt * w_inc, w_star)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    plt.figure(figsize=(7, 3.5))
    plt.plot(t, v, '-k', linewidth=2)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('$v$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
