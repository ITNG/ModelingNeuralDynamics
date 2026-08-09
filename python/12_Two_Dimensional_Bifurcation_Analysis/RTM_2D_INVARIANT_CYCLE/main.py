import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

c = 1.
g_k = 80.
g_na = 100.
g_l = 0.1
v_k = -100.
v_na = 50.
v_l = -67.

xlim, ylim = [-75, -50], [0, 0.15]
frame_xlim, frame_ylim = [-105, 55], [-0.05, 0.8]


def alpha_m(v):
    return 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))


def beta_m(v):
    return 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)


def alpha_n(v):
    return 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))


def beta_n(v):
    return 0.5 * exp(-(v + 57) / 40)


def alpha_h(v):
    return 0.128 * exp(-(v + 50) / 18)


def beta_h(v):
    return 4. / (1 + exp(-(v + 27) / 5))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def simulate(v0, n0, i_ext, t_final, dt=0.005):
    '''reduced RTM model: m=m_inf(v), h=1-n'''
    m_steps = round(t_final / dt)
    v, n = np.zeros(m_steps + 1), np.zeros(m_steps + 1)
    v[0], n[0] = v0, n0
    for k in range(m_steps):
        h = 1 - n[k]
        v_inc = (g_k * n[k] ** 4 * (v_k - v[k]) + g_na * m_inf(v[k]) ** 3 * h * (v_na - v[k])
                 + g_l * (v_l - v[k]) + i_ext) / c
        n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
        v_tmp = v[k] + dt / 2 * v_inc
        n_tmp = n[k] + dt / 2 * n_inc
        h_tmp = 1 - n_tmp
        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_inf(v_tmp) ** 3 * h_tmp * (v_na - v_tmp)
                 + g_l * (v_l - v_tmp) + i_ext) / c
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        v[k + 1] = v[k] + dt * v_inc
        n[k + 1] = n[k] + dt * n_inc
    return v, n


def find_root(f, a, b, tol):
    while b - a > tol:
        c_ = (a + b) / 2
        if f(c_) * f(a) > 0:
            a = c_
        else:
            b = c_
    return c_


def mark_crossings(ax, v, n, condition, ax_xlim, ax_ylim, color='k'):
    for k in np.where(condition)[0]:
        vec = [v[k + 1] - v[k], n[k + 1] - n[k]]
        draw_arrow(ax, ax_xlim, ax_ylim, v[k], n[k], vec, epsilon=0.1, width=2, color=color)


if __name__ == "__main__":

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))

    # panel 1: i_ext=1 -- limit cycle, spiraling in from a nearby point
    v, n = simulate(-75.582445796204553, 0.005935897840905, i_ext=1., t_final=30.)
    mark_crossings(axes[0], v, n, (v[:-1] < 0) & (v[1:] > 0), frame_xlim, frame_ylim)
    mark_crossings(axes[0], v, n, (v[:-1] > -30) & (v[1:] < -30), frame_xlim, frame_ylim)
    axes[0].plot(v, n, color='k', linewidth=2)
    axes[0].set_xlim(*frame_xlim)
    axes[0].set_ylim(*frame_ylim)
    axes[0].set_box_aspect(1)
    axes[0].set_xlabel('$v$ [mV]')
    axes[0].set_ylabel('$n$')

    # i_ext=0 -- the two fixed points that bound the invariant cycle
    i_ext = 0.

    def f(v):
        # note: h_inf(v) here, not the 1-n reduction simulate() uses --
        # matches matlab/12/RTM_2D_INVARIANT_CYCLE/make_figure.m exactly
        return g_k * n_inf(v) ** 4 * (v_k - v) + g_na * m_inf(v) ** 3 * h_inf(v) * (v_na - v) \
            + g_l * (v_l - v) + i_ext

    v_star = find_root(f, -63, -62, 1e-14)
    n_star = n_inf(v_star)
    v_0 = find_root(f, -75, -65, 1e-12)
    n_0 = n_inf(v_0)

    # panel 2: spiral converging onto the invariant cycle from just outside
    v, n = simulate(v_star + 0.5, n_star + 0.005, i_ext=0., t_final=300.)
    mark_crossings(axes[1], v, n, (v[:-1] < 0) & (v[1:] > 0), frame_xlim, frame_ylim)
    mark_crossings(axes[1], v, n, (v[:-1] > -30) & (v[1:] < -30), frame_xlim, frame_ylim)
    axes[1].plot(v, n, color='k', linewidth=2)
    axes[1].set_xlim(*frame_xlim)
    axes[1].set_ylim(*frame_ylim)
    axes[1].set_box_aspect(1)
    axes[1].set_xlabel('$v$ [mV]')

    # panel 3: zoomed view of the same cycle, plus one spiraling out from
    # just inside v_star (the unstable fixed point) toward the cycle
    axes[2].plot(v, n, color='k', linewidth=2)
    ind = np.where((v[:-1] < -55) & (v[1:] >= -55))[0][0]
    vec = [v[ind + 1] - v[ind], n[ind + 1] - n[ind]]
    draw_arrow(axes[2], xlim, ylim, v[ind], n[ind], vec, epsilon=0.1, width=2)
    ind = np.where((v[:-1] < -70) & (v[1:] >= -70))[0][0]
    vec = [v[ind + 1] - v[ind], n[ind + 1] - n[ind]]
    draw_arrow(axes[2], xlim, ylim, v[ind], n[ind], vec, epsilon=0.1, width=2)

    v2, n2 = simulate(v_star - 0.5, n_star - 0.005, i_ext=0., t_final=300.)
    mark_crossings(axes[2], v2, n2, (v2[:-1] > -65) & (v2[1:] < -65), xlim, ylim, color='r')
    axes[2].plot(v2, n2, color='r', linewidth=2)
    axes[2].plot(v_star, n_star, 'ok', markersize=8, markerfacecolor='w')
    axes[2].plot(v_0, n_0, 'ok', markersize=8, markerfacecolor='k')
    axes[2].set_xlim(*xlim)
    axes[2].set_ylim(*ylim)
    axes[2].set_box_aspect(1)
    axes[2].set_xlabel('$v$ [mV]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
