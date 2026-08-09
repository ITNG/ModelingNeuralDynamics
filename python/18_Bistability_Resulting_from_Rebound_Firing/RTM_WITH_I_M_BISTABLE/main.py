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


def simulate(v0, m0, h0, n0, w0):
    '''Heun/RK2 integration (plain floats) of the RTM model with an added
    M-current, m quasi-static'''
    v, m, h, n, w = v0, m0, h0, n0, w0
    v_trace = np.empty(m_steps + 1)
    v_trace[0] = v
    for k in range(m_steps):
        v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v)
                  + g_m * w * (v_k - v) + i_ext) / c
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
        w_inc = (w_inf(v) - w) / tau_w(v)

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        w_tmp = w + dt05 * w_inc

        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_l * (v_l - v_tmp) + g_m * w_tmp * (v_k - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        # matlab's second-stage w_inc divides by w(k), not w_tmp -- kept
        # as-is for a faithful port (see ch17/RTM_WITH_M_CURRENT_F_I)
        w_inc = (w_inf(v_tmp) - w) / tau_w(v_tmp)

        v = v + dt * v_inc
        m = m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc
        w = w + dt * w_inc
        v_trace[k + 1] = v
    return v_trace


v_rest = simulate(v_star, m_star, h_star, n_star, w_star)
v_fire = simulate(v_star + 1, m_star, h_star, n_star, w_star)

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
