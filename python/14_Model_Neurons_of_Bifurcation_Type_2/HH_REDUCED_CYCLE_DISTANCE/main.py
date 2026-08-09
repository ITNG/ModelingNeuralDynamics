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


def alpha_m(v):
    # removable singularity at v=-45 (0/0) -- matlab's alpha_m.m branches
    # to the L'Hopital limit (1) there instead of evaluating it
    if abs(v + 45) > 1e-8:
        return (v + 45) / 10. / (1 - exp(-(v + 45) / 10))
    return 1.0


def alpha_n(v):
    return 0.01 * (-60. - v) / (exp((-60 - v) / 10) - 1)


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def f(v):
    '''dv/dt of the reduced (m=m_inf(v), h=0.83-n_inf(v)) HH model, i_ext=0'''
    return g_na * m_inf(v) ** 3 * (0.83 - n_inf(v)) * (v_na - v) \
        + g_k * n_inf(v) ** 4 * (v_k - v) + g_l * (v_l - v)


def find_fixed_point(i_ext):
    '''bisect f(v)+i_ext=0'''
    v_left, v_right = -100., 50.
    while v_right - v_left > 1e-10:
        v_c = (v_left + v_right) / 2
        if (f(v_c) + i_ext) * (f(v_left) + i_ext) > 0:
            v_left = v_c
        else:
            v_right = v_c
    return (v_left + v_right) / 2


def simulate(v0, n0, i_ext, t_final, dt, direction):
    '''Heun/RK2 integration of the reduced HH model.

    direction=+1 integrates forward in time, direction=-1 backward
    (used to trace the unstable manifold near the repelling cycle).
    '''
    dt05 = dt / 2
    m_steps = round(t_final / dt)
    v = np.zeros(m_steps + 1)
    n = np.zeros(m_steps + 1)
    v[0], n[0] = v0, n0

    for k in range(m_steps):
        m_k = m_inf(v[k])
        h_k = 0.83 - n[k]
        v_inc = (g_na * m_k ** 3 * h_k * (v_na - v[k]) + g_k * n[k] ** 4 * (v_k - v[k])
                  + g_l * (v_l - v[k]) + i_ext) / c
        n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]

        v_tmp = v[k] + direction * dt05 * v_inc
        n_tmp = n[k] + direction * dt05 * n_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = 0.83 - n_tmp

        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v[k + 1] = v[k] + direction * dt * v_inc
        n[k + 1] = n[k] + direction * dt * n_inc

    return v, n


t_final = 1000.
dt = 0.01

i_ext_vec = [5.5, 5.4, 5.3]
panels = []
for i_ext in i_ext_vec:
    v_c = find_fixed_point(i_ext)
    n_c = n_inf(v_c)

    v_attr, n_attr = simulate(20., n_c, i_ext, t_final, dt, direction=1)
    tail = round((t_final - 20) / dt)
    v_attr, n_attr = v_attr[tail:], n_attr[tail:]

    v_rep, n_rep = simulate(v_c + 0.001, n_c, i_ext, t_final, dt, direction=-1)
    tail = round((t_final - 15) / dt)
    v_rep, n_rep = v_rep[tail:], n_rep[tail:]

    panels.append((i_ext, v_c, n_c, v_attr, n_attr, v_rep, n_rep))

if __name__ == "__main__":

    fig, ax = plt.subplots(1, 3, figsize=(12, 4.5))

    for a, (i_ext, v_c, n_c, v_attr, n_attr, v_rep, n_rep) in zip(ax, panels):
        a.plot(v_attr, n_attr, '-k', linewidth=2)
        a.plot(v_c, n_c, '.', markersize=25)
        a.plot(v_rep, n_rep, ':r', linewidth=2)
        a.set_xlim(-67.5, -65.5)
        a.set_ylim(0.3575, 0.3675)
        a.set_xticks([-67, -66])
        a.set_yticks([0.36, 0.365])
        a.set_xlabel('$v$')
        a.set_title(f'$I={i_ext}$')
        a.set_box_aspect(1)

    ax[0].set_ylabel('$n$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
