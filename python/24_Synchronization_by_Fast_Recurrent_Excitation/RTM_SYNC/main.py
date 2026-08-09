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

i_ext = 0.3

N = 30  # number of neurons in the network

t_final = 200.
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


v = -70. * np.ones(N)
m = m_inf(v)
h = h_inf(v)
n = n_inf(v)

t_spikes = []
i_spikes = []

for k in range(1, m_steps + 1):
    v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v) + i_ext) / c
    n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
    h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h

    v_tmp = v + dt05 * v_inc
    h_tmp = h + dt05 * h_inc
    n_tmp = n + dt05 * n_inc
    m_tmp = m_inf(v_tmp)

    v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
              + g_l * (v_l - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

    v_old = v
    v = v + dt * v_inc
    h = h + dt * h_inc
    n = n + dt * n_inc
    m = m_inf(v)

    ind = np.where((v_old >= -20) & (v < -20))[0]
    if len(ind) > 0:
        t_spikes.extend(((k - 1) * dt * (-20 - v[ind]) + k * dt * (v_old[ind] + 20))
                         / (v_old[ind] - v[ind]))
        i_spikes.extend(ind)

t_spikes = np.array(t_spikes)
i_spikes = np.array(i_spikes)

if __name__ == "__main__":

    plt.figure(figsize=(7, 4))
    plt.plot(t_spikes, i_spikes + 1, '.r', markersize=10)
    plt.xlim(0, t_final)
    plt.ylim(0, N + 1)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('neuron #')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
