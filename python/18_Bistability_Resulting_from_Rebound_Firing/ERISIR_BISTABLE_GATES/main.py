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

i_ext = 6.9

t_final = 40.
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


def f(v):
    return (g_k * n_inf(v) ** 2 * (v_k - v) + g_na * m_inf(v) ** 3 * h_inf(v) * (v_na - v)
            + g_l * (v_l - v) + i_ext) / c


v_vec = -100 + np.arange(500001) / 500000 * 150
with np.errstate(invalid='ignore'):
    f_vec = f(v_vec)
ind = np.where(f_vec[:-1] * f_vec[1:] < 0)[0].min()
v_star = (v_vec[ind] + v_vec[ind + 1]) / 2
m_star = m_inf(v_star)
h_star = h_inf(v_star)
n_star = n_inf(v_star)

v = np.empty(m_steps + 1)
m = np.empty(m_steps + 1)
h = np.empty(m_steps + 1)
n = np.empty(m_steps + 1)
v[0], m[0], h[0], n[0] = v_star + 3, m_star, h_star, n_star

for k in range(m_steps):
    v_inc = (g_k * n[k] ** 2 * (v_k - v[k]) + g_na * m[k] ** 3 * h[k] * (v_na - v[k]) + g_l * (v_l - v[k]) + i_ext) / c
    n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
    h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]

    v_tmp = v[k] + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h[k] + dt05 * h_inc
    n_tmp = n[k] + dt05 * n_inc

    v_inc = (g_k * n_tmp ** 2 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
              + g_l * (v_l - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

    v[k + 1] = v[k] + dt * v_inc
    m[k + 1] = m_inf(v[k + 1])
    h[k + 1] = h[k] + dt * h_inc
    n[k + 1] = n[k] + dt * n_inc

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    fig, ax = plt.subplots(3, figsize=(7, 8), sharex=True)

    ax[0].plot(t, m_inf(v), '--b', linewidth=2)
    ax[0].plot([0, t_final], [m_star, m_star], '-r', linewidth=2)
    ax[0].set_ylabel('$m$')
    ax[0].set_ylim(0, 1)

    ax[1].plot(t, h, '-k', linewidth=2)
    ax[1].plot(t, h_inf(v), '--b', linewidth=2)
    ax[1].plot([0, t_final], [h_star, h_star], '-r', linewidth=2)
    ax[1].set_ylabel('$h$')
    ax[1].set_ylim(0, 1)

    ax[2].plot(t, n, '-k', linewidth=2)
    ax[2].plot(t, n_inf(v), '--b', linewidth=2)
    ax[2].plot([0, t_final], [n_star, n_star], '-r', linewidth=2)
    ax[2].set_ylabel('$n$')
    ax[2].set_xlabel('$t$ [ms]')
    ax[2].set_ylim(0, 1)
    ax[2].set_xlim(0, t_final)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
