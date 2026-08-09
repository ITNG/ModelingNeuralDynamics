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
    return out if abs(v + 45) > 1e-8 else 1.0


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
    '''dv/dt at the given i_ext, for the resting-state bisection'''
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


def simulate(v0, m0, h0, n0):
    '''Heun/RK2 integration (plain floats) of the full (v,m,h,n) HH model'''
    v, m, h, n = v0, m0, h0, n0
    v_trace = np.empty(m_steps + 1)
    v_trace[0] = v
    for k in range(m_steps):
        v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext) / c
        m_inc = alpha_m(v) * (1 - m) - beta_m(v) * m
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n

        v_tmp = v + dt05 * v_inc
        m_tmp = m + dt05 * m_inc
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc

        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        m_inc = alpha_m(v_tmp) * (1 - m_tmp) - beta_m(v_tmp) * m_tmp
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v = v + dt * v_inc
        m = m + dt * m_inc
        h = h + dt * h_inc
        n = n + dt * n_inc
        v_trace[k + 1] = v
    return v_trace


v_rest = simulate(v_star, m_star, h_star, n_star)
v_fire = simulate(v_star + 5, m_star, h_star, n_star)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    fig, ax = plt.subplots(2, figsize=(7, 6), sharex=True)
    ax[0].plot(t, v_rest, '-k', linewidth=2)
    ax[0].set_ylabel('$v$ [mV]')
    ax[0].set_xlim(0, t_final)
    ax[0].set_ylim(-100, 50)

    ax[1].plot(t, v_fire, '-k', linewidth=2)
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$v$ [mV]')
    ax[1].set_xlim(0, t_final)
    ax[1].set_ylim(-100, 50)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
