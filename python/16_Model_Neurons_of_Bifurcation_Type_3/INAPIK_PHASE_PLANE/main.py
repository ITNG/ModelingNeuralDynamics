import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

c = 1.
g_na = 20.
g_k = 10.
g_l = 8.
v_na = 60.
v_k = -90.
v_l = -80.
tau_n = 0.15


def m_inf(v):
    return 1. / (1 + exp((-20 - v) / 15))


def m_inf_p(v):
    return -1. / (1 + exp((-20 - v) / 15)) ** 2 * exp((-20 - v) / 15) * (-1 / 15)


def n_inf(v):
    return 1. / (1 + exp((-25 - v) / 5))


def n_inf_p(v):
    return -1. / (1 + exp((-25 - v) / 5)) ** 2 * exp((-25 - v) / 5) * (-1 / 5)


def f(v):
    '''zeros of this function (offset by I) are the fixed points'''
    return g_na * m_inf(v) * (v_na - v) + g_k * n_inf(v) * (v_k - v) + g_l * (v_l - v)


w_grid = -100 + np.arange(10001) / 10000 * 150
f_grid = f(w_grid)


def find_fixed_points(I):
    fi = f_grid + I
    roots = []
    for j in np.where(fi[:-1] * fi[1:] <= 0)[0]:
        w_low, w_high = w_grid[j], w_grid[j + 1]
        while w_high - w_low > 1e-12:
            w_c = (w_low + w_high) / 2
            if (f(w_c) + I) * (f(w_high) + I) <= 0:
                w_low = w_c
            else:
                w_high = w_c
        roots.append((w_low + w_high) / 2)
    return roots


def jacobian(v_c, n_c):
    j00 = g_na * m_inf_p(v_c) * (v_na - v_c) - g_na * m_inf(v_c) - g_k * n_c - g_l
    j01 = g_k * (v_k - v_c)
    j10 = n_inf_p(v_c) / tau_n
    j11 = -1 / tau_n
    return np.array([[j00, j01], [j10, j11]])


def classify_fixed_points(I):
    out = {'black': [], 'blue': [], 'red': [], 'green': [], 'magenta': []}
    for v_c in find_fixed_points(I):
        n_c = n_inf(v_c)
        e = np.linalg.eigvals(jacobian(v_c, n_c))
        if abs(e[0].imag) < 1e-12 and e[0].real < 0 and e[1].real < 0:
            out['black'].append((v_c, n_c))
        if abs(e[0].imag) < 1e-12 and e[0].real > 0 and e[1].real > 0:
            out['blue'].append((v_c, n_c))
        if abs(e[0].imag) > 1e-12 and e[0].real < 0:
            out['red'].append((v_c, n_c))
        if abs(e[0].imag) > 1e-12 and e[0].real > 0:
            out['green'].append((v_c, n_c))
        if abs(e[0].imag) < 1e-12 and e[0].real * e[1].real < 0:
            out['magenta'].append((v_c, n_c))
    return out


def simulate(v0, n0, i_ext, t_final, dt):
    dt05 = dt / 2
    m_steps = round(t_final / dt)
    v = np.zeros(m_steps + 1)
    n = np.zeros(m_steps + 1)
    v[0], n[0] = v0, n0

    for k in range(m_steps):
        m = m_inf(v[k])
        v_inc = (g_na * m * (v_na - v[k]) + g_k * n[k] * (v_k - v[k]) + g_l * (v_l - v[k]) + i_ext) / c
        n_inc = (n_inf(v[k]) - n[k]) / tau_n

        v_tmp = v[k] + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        n_tmp = n[k] + dt05 * n_inc

        v_inc = (g_na * m_tmp * (v_na - v_tmp) + g_k * n_tmp * (v_k - v_tmp) + g_l * (v_l - v_tmp) + i_ext) / c
        n_inc = (n_inf(v_tmp) - n_tmp) / tau_n

        v[k + 1] = v[k] + dt * v_inc
        n[k + 1] = n[k] + dt * n_inc

    return v, n


def _plot_fixed_points(ax, fp):
    for v_c, n_c in fp['black']:
        ax.plot(v_c, n_c, '.k', markersize=15)
    for v_c, n_c in fp['blue']:
        ax.plot(v_c, n_c, '.b', markersize=15)
    for v_c, n_c in fp['red']:
        ax.plot(v_c, n_c, '.r', markersize=15)
    for v_c, n_c in fp['green']:
        ax.plot(v_c, n_c, 'og', markersize=6, markerfacecolor='none')
    for v_c, n_c in fp['magenta']:
        ax.plot(v_c, n_c, 'om', markersize=6, markerfacecolor='none')


A, B, C, D = -75, 0, -0.15, 0.85


def _arrow(ax, v, n, t0, dt, epsilon=0.07, theta=0.):
    '''arrow position/direction from the *full* trajectory at time t0 --
    for the limit-cycle panels t0 can fall before the displayed tail
    window, landing on an earlier pass through the same periodic orbit'''
    k0 = round(t0 / dt)
    x0, y0 = v[k0], n[k0]
    vec = np.array([(v[k0 + 1] - v[k0 - 1]) / (2 * dt),
                     (n[k0 + 1] - n[k0 - 1]) / (2 * dt)])
    if theta:
        rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        vec = rot @ vec
    draw_arrow(ax, (A, B), (C, D), x0, y0, vec, epsilon=epsilon, width=2, color='k')


panels = []

# I=-5: settles at a stable fixed point
v1, n1 = simulate(-50., 0., i_ext=-5., t_final=2.5, dt=0.001)
panels.append(dict(v_full=v1, n_full=n1, v=v1, n=n1, dt=0.001, i_ext=-5., title=r'$I=-5$',
                    arrows=[(0.4, 0.), (0.7, 0.), (1.8, 0.002)]))

# I=-1.4: settles at a stable fixed point (near the SNIC)
v2, n2 = simulate(-52.5, 0., i_ext=-1.4, t_final=3., dt=0.001)
panels.append(dict(v_full=v2, n_full=n2, v=v2, n=n2, dt=0.001, i_ext=-1.4, title=r'$I=-1.4$',
                    arrows=[(0.5, 0.), (0.9, 0.)]))

# I=2: settles onto the stable limit cycle -- only the last 1.2ms is plotted
t_final = 20.
dt = 0.001
v3, n3 = simulate(-30., 0.2, i_ext=2., t_final=t_final, dt=dt)
tail = round((t_final - 1.2) / dt)
panels.append(dict(v_full=v3, n_full=n3, v=v3[tail:], n=n3[tail:], dt=dt, i_ext=2., title=r'$I=2$',
                    arrows=[(t_final - 1.6, 0.), (t_final - 1.9, 0.)]))

# I=4.4: settles onto the stable limit cycle -- only the last 2ms plotted
v4, n4 = simulate(-30., 0.2, i_ext=4.4, t_final=t_final, dt=dt)
tail = round((t_final - 2) / dt)
panels.append(dict(v_full=v4, n_full=n4, v=v4[tail:], n=n4[tail:], dt=dt, i_ext=4.4, title=r'$I=4.4$',
                    arrows=[(t_final - 1, 0.), (t_final - 1.3, 0.)]))

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 2, figsize=(9, 9))

    for ax, p in zip(axes.flat, panels):
        ax.plot(p['v'], p['n'], '-k', linewidth=2)
        fp = classify_fixed_points(p['i_ext'])
        _plot_fixed_points(ax, fp)
        for t0, theta in p['arrows']:
            _arrow(ax, p['v_full'], p['n_full'], t0, p['dt'], theta=theta)
        ax.set_xlim(A, B)
        ax.set_ylim(C, D)
        ax.set_box_aspect(1)
        ax.set_title(p['title'])

    axes[1, 0].set_xlabel('$v$')
    axes[1, 0].set_ylabel('$n$')
    axes[1, 1].set_xlabel('$v$')
    axes[0, 0].set_ylabel('$n$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
