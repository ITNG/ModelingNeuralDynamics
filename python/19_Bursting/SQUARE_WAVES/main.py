import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_na = 20.
g_k = 10.
g_l = 8.
v_na = 60.
v_k = -90.
v_l = -80.
tau_n = 0.15

i_ext = 7.

g_k_slow = 5.
tau_n_slow = 20.

t_final = 100.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


def m_inf(v):
    return 1. / (1 + exp((-20 - v) / 15))


def n_inf(v):
    return 1. / (1 + exp((-25 - v) / 5))


def n_slow_inf(v):
    return 1. / (1 + exp((-20 - v) / 5))


v = np.zeros(m_steps + 1)
m = np.zeros(m_steps + 1)
n = np.zeros(m_steps + 1)
n_slow = np.zeros(m_steps + 1)
v[0], m[0], n[0], n_slow[0] = -70., m_inf(-70.), 0.6, 0.

for k in range(m_steps):
    v_inc = (g_na * m[k] * (v_na - v[k]) + g_k * n[k] * (v_k - v[k])
              + g_k_slow * n_slow[k] * (v_k - v[k]) + g_l * (v_l - v[k]) + i_ext) / c
    n_inc = (n_inf(v[k]) - n[k]) / tau_n
    n_slow_inc = (n_slow_inf(v[k]) - n_slow[k]) / tau_n_slow

    v_tmp = v[k] + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    n_tmp = n[k] + dt05 * n_inc
    n_slow_tmp = n_slow[k] + dt05 * n_slow_inc

    v_inc = (g_na * m_tmp * (v_na - v_tmp) + g_k * n_tmp * (v_k - v_tmp)
              + g_k_slow * n_slow_tmp * (v_k - v_tmp) + g_l * (v_l - v_tmp) + i_ext) / c
    n_inc = (n_inf(v_tmp) - n_tmp) / tau_n
    n_slow_inc = (n_slow_inf(v_tmp) - n_slow_tmp) / tau_n_slow

    v[k + 1] = v[k] + dt * v_inc
    m[k + 1] = m_inf(v[k + 1])
    n[k + 1] = n[k] + dt * n_inc
    n_slow[k + 1] = n_slow[k] + dt * n_slow_inc

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    plt.figure(figsize=(7, 3.5))
    plt.plot(t, v, '-k', linewidth=1)

    # highlight the slow (quasi-steady) portions of the trajectory, where
    # the central difference dv/dt is small
    v_left, v_right = v[:-2], v[2:]
    ind = np.where(np.abs(v_right - v_left) / dt < 1)[0]
    plt.plot(t[ind + 1], v[ind + 1], '-r', linewidth=2)

    plt.xlabel('$t$ [ms]')
    plt.ylabel('$v$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
