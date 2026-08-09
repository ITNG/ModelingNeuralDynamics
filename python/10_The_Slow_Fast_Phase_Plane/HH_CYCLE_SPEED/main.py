import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.0
g_k = 36.0
g_na = 120.0
g_l = 0.3
v_k = -82.0
v_na = 45.0
v_l = -59.0
i_ext = 10.0


def alpha_m(v):
    # removable singularity at v=-45 (0/0) -- matlab/10/HH_CYCLE_SPEED/alpha_m.m
    # branches to the L'Hopital limit (1) there instead of evaluating it
    if abs(v + 45) > 1e-8:
        return (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))
    return 1.0


def alpha_n(v):
    return 0.01 * (-60.0 - v) / (exp((-60 - v) / 10) - 1)


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


# v-nullcline: for each v, the n where the total membrane current is zero
# (found by bisection, since it's not solvable for n in closed form)
v_vec = np.arange(-100, 51)
n_vec = np.zeros_like(v_vec, dtype=float)
for ij, v in enumerate(v_vec):
    n_l, n_r = 0.0, 1.0
    while n_r - n_l > 1e-10:
        n = (n_l + n_r) / 2
        i_ion = g_na * m_inf(v) ** 3 * (0.83 - n) * (v_na - v) \
            + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext
        if i_ion < 0:
            n_r = n
        else:
            n_l = n
    n_vec[ij] = (n_l + n_r) / 2
    if n_vec[ij] > 1 - 1e-10:
        n_vec[ij] = 2
    if n_vec[ij] < 1e-10:
        n_vec[ij] = -1

t_final = 150.
dt = 0.001
dt05 = dt / 2
m_steps = round(t_final / dt)

v = np.zeros(m_steps + 1)
n = np.zeros(m_steps + 1)
speed = np.zeros(m_steps)

v[0] = -70.
n[0] = 0.

for k in range(m_steps):
    m_k = m_inf(v[k])
    h_k = 0.83 - n[k]
    v_inc = (g_na * m_k ** 3 * h_k * (v_na - v[k]) + g_k * n[k] ** 4 * (v_k - v[k])
              + g_l * (v_l - v[k]) + i_ext) / c
    n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
    v_tmp = v[k] + dt05 * v_inc
    n_tmp = n[k] + dt05 * n_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = 0.83 - n_tmp
    v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
              + g_l * (v_l - v_tmp) + i_ext) / c
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
    v[k + 1] = v[k] + dt * v_inc
    n[k + 1] = n[k] + dt * n_inc
    speed[k] = np.sqrt((v_inc / 150) ** 2 + (n_inc / 0.35) ** 2)

n_start = round(2 / 3 * m_steps)
v = v[n_start:m_steps]
n = n[n_start:m_steps]
speed = speed[n_start:m_steps]

if __name__ == "__main__":

    plt.figure(figsize=(7, 7))
    # offset by 0.001, same as matlab -- alpha_n has its own (unguarded)
    # removable singularity at v=-60, which is otherwise hit exactly
    v_line = v_vec.astype(float) + 0.001
    plt.plot(v_line, n_inf(v_line), color='k', linewidth=3)
    plt.plot(v_vec, n_vec, color='r', linewidth=3)

    plt.plot(v, n, color='b', linewidth=2)
    fast = speed > speed.max() * 0.02
    plt.plot(v[fast], n[fast], '.b', markersize=10)
    plt.plot(v[~fast], n[~fast], '.g', markersize=10)

    plt.xlim(-100, 50)
    plt.ylim(0, 1)
    plt.xlabel('$v$ [mV]')
    plt.ylabel('$n$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
