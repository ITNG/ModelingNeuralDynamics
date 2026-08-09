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

g_h = 1.
v_h = -32.9

i_ext = -3.19

t_final = 3000.
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


def r_inf(v):
    return 1. / (1 + exp((v + 84) / 10.2))


def tau_r(v):
    return 1. / (exp(-14.59 - 0.086 * v) + exp(-1.87 + 0.0701 * v))


def f(v):
    return g_na * m_inf(v) ** 3 * h_inf(v) * (v_na - v) + g_k * n_inf(v) ** 4 * (v_k - v) \
        + g_l * (v_l - v) + g_h * r_inf(v) * (v_h - v) + i_ext


v_grid = -100 + np.arange(100001) / 100000 * 150
f_grid = f(v_grid)
ind = np.where(f_grid[:-1] * f_grid[1:] <= 0)[0]
v_left = v_grid[ind].min()
v_right = v_grid[ind + 1].min()
while v_right - v_left > 1e-12:
    v_c = (v_left + v_right) / 2
    if f(v_c) * f(v_left) <= 0:
        v_right = v_c
    else:
        v_left = v_c
v_star = (v_left + v_right) / 2
m_star = m_inf(v_star)
h_star = h_inf(v_star)
n_star = n_inf(v_star)
r_star = r_inf(v_star)


def simulate(v0, m0, h0, n0, r0):
    '''Heun/RK2 integration (plain floats) of the RTM model with an added
    I_h current, m quasi-static'''
    v, m, h, n, r = v0, m0, h0, n0, r0
    v_trace = np.empty(m_steps + 1)
    v_trace[0] = v
    for k in range(m_steps):
        v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v)
                  + g_h * r * (v_h - v) + i_ext) / c
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
        r_inc = (r_inf(v) - r) / tau_r(v)

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        r_tmp = r + dt05 * r_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + g_h * r_tmp * (v_h - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        r_inc = (r_inf(v_tmp) - r_tmp) / tau_r(v_tmp)

        v = v + dt * v_inc
        m = m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc
        r = r + dt * r_inc
        v_trace[k + 1] = v
    return v_trace


v_rest = simulate(v_star + 0.01, m_star, h_star, n_star, r_star)
v_fire = simulate(v_star + 1, m_star, h_star, n_star, r_star)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    fig, ax = plt.subplots(2, figsize=(7, 6), sharex=True)
    ax[0].plot(t, v_rest, '-k', linewidth=2)
    ax[0].set_ylabel('$v$ [mV]')

    ax[1].plot(t, v_fire, '-k', linewidth=2)
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$v$ [mV]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
