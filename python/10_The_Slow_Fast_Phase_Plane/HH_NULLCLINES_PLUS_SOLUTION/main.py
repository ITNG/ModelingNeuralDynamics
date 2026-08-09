from scipy.integrate import odeint
import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

c = 1.0
g_k = 36.0
g_na = 120.0
g_l = 0.3
v_k = -82.0
v_na = 45.0
v_l = -59.0
i_ext = 10.0


def alpha_m(v):
    # removable singularity at v=-45 (0/0) -- matlab's alpha_m.m branches
    # to the L'Hopital limit (1) there instead of evaluating it
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


def derivative(x0, t):
    v, n = x0
    m = m_inf(v)
    h = 0.83 - n
    dv = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v)
          + g_l * (v_l - v) + i_ext) / c
    dn = alpha_n(v) * (1 - n) - beta_n(v) * n
    return [dv, dn]


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

t_final = 50.
dt = 0.01
x0 = [-70., 0.]

t = np.arange(0, t_final, dt)
sol = odeint(derivative, x0, t)
v, n = sol[:, 0], sol[:, 1]

if __name__ == "__main__":

    fig, ax = plt.subplots(figsize=(7, 7))

    # offset by 0.001, same as matlab -- alpha_n has its own (unguarded)
    # removable singularity at v=-60, which is otherwise hit exactly
    v_line = v_vec.astype(float) + 0.001
    ax.plot(v_line, n_inf(v_line), color='k', linewidth=3)
    ax.plot(v_vec, n_vec, color='r', linewidth=3)
    ax.plot(v, n, color='b', linewidth=2)

    for t_arrow in [0.57, 2.5, 6.5, 12.75, 13.2]:
        i = round(t_arrow / dt)
        vec = sol[i] - sol[i - 11]
        draw_arrow(ax, [-100, 50], [0, 1], v[i], n[i], vec, epsilon=0.05, color='b')

    ax.set_xlim(-100, 50)
    ax.set_ylim(0, 1)
    ax.set_xlabel('$v$ [mV]')
    ax.set_ylabel('$n$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
