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
f_vec = f(v_vec)
ind = np.where(f_vec[:-1] * f_vec[1:] < 0)[0].min()
v_star = (v_vec[ind] + v_vec[ind + 1]) / 2
m_star = m_inf(v_star)
h_star = h_inf(v_star)
n_star = n_inf(v_star)


def simulate(v0, m0, h0, n0):
    '''Heun/RK2 integration (plain floats) of the Erisir model, m
    quasi-static, h and n dynamic'''
    v, m, h, n = v0, m0, h0, n0
    v_trace = np.empty(m_steps + 1)
    v_trace[0] = v
    for k in range(m_steps):
        v_inc = (g_k * n ** 2 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v) + i_ext) / c
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc

        v_inc = (g_k * n_tmp ** 2 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v = v + dt * v_inc
        m = m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc
        v_trace[k + 1] = v
    return v_trace


v_rest = simulate(v_star, m_star, h_star, n_star)
v_fire = simulate(v_star + 3, m_star, h_star, n_star)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    fig, ax = plt.subplots(2, figsize=(7, 6), sharex=True)
    ax[0].plot(t, v_rest, '-k', linewidth=2)
    ax[0].set_ylabel('$v$ [mV]')
    ax[0].set_xlim(0, t_final)
    ax[0].set_ylim(-95, 55)

    ax[1].plot(t, v_fire, '-k', linewidth=2)
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$v$ [mV]')
    ax[1].set_xlim(0, t_final)
    ax[1].set_ylim(-95, 55)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
