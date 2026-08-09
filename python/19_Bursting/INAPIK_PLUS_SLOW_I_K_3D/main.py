import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_na = 20.
g_k = 10.
g_k_slow = 5.
g_l = 8.
v_na = 60.
v_k = -90.
v_l = -80.
tau_n = 0.15
tau_n_slow = 20.

i_ext = 7.

t_final = 1500.
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

k_vec = []
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

    # a full pass around the limit cycle is complete each time n_slow
    # drops back below 0.02
    if n_slow[k] > 0.02 and n_slow[k + 1] <= 0.02:
        k_vec.append(k)

# only the last full passage through the limit cycle, so the plot shows
# the cycle itself but not the transient approach to it
k1, k2 = k_vec[-2], k_vec[-1]
v_cyc = v[k1:k2 + 1]
n_cyc = n[k1:k2 + 1]
n_slow_cyc = n_slow[k1:k2 + 1]

if __name__ == "__main__":

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection='3d')
    ax.plot(v_cyc, n_cyc, n_slow_cyc, '-k', linewidth=2)
    ax.set_xlabel('$v$')
    ax.set_ylabel('$n$')
    ax.set_zlabel(r'$n_{\rm slow}$')
    ax.set_xlim(-70, 0)
    ax.set_ylim(0, 0.7)
    ax.set_zlim(0, 0.06)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
