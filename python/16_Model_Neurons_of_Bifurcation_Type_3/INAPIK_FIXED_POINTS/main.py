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


v_grid = -100 + np.arange(3001) / 3000 * 150
f_grid = f(v_grid)


def find_fixed_points(I):
    '''bracket sign changes of f+I on v_grid, then bisect each one'''
    fi = f_grid + I
    roots = []
    for j in np.where(fi[:-1] * fi[1:] <= 0)[0]:
        v_low, v_high = v_grid[j], v_grid[j + 1]
        while v_high - v_low > 1e-12:
            v_c = (v_low + v_high) / 2
            if (f(v_c) + I) * (f(v_high) + I) <= 0:
                v_low = v_c
            else:
                v_high = v_c
        roots.append((v_low + v_high) / 2)
    return roots


def jacobian(v_c, n_c):
    j00 = g_na * m_inf_p(v_c) * (v_na - v_c) - g_na * m_inf(v_c) - g_k * n_c - g_l
    j01 = g_k * (v_k - v_c)
    j10 = n_inf_p(v_c) / tau_n
    j11 = -1 / tau_n
    return np.array([[j00, j01], [j10, j11]])


i_ext_vec = -4 + np.arange(1001) / 1000 * 12

points = {'black': [], 'blue': [], 'red': [], 'green': [], 'magenta': []}
for I in i_ext_vec:
    for v_c in find_fixed_points(I):
        n_c = n_inf(v_c)
        e = np.linalg.eigvals(jacobian(v_c, n_c))
        if abs(e[0].imag) < 1e-12 and e[0].real < 0 and e[1].real < 0:
            points['black'].append((I, v_c))
        if abs(e[0].imag) < 1e-12 and e[0].real > 0 and e[1].real > 0:
            points['blue'].append((I, v_c))
        if abs(e[0].imag) > 1e-12 and e[0].real < 0:
            points['red'].append((I, v_c))
        if abs(e[0].imag) > 1e-12 and e[0].real > 0:
            points['green'].append((I, v_c))
        if abs(e[0].imag) < 1e-12 and e[0].real * e[1].real < 0:
            points['magenta'].append((I, v_c))

if __name__ == "__main__":

    all_v = [v for pts in points.values() for _, v in pts]
    maxv_c, minv_c = max(all_v + [-100]), min(all_v + [50])

    plt.figure(figsize=(7, 7))
    for color, key in [('k', 'black'), ('b', 'blue'), ('r', 'red')]:
        if points[key]:
            i_pts, v_pts = zip(*points[key])
            plt.plot(i_pts, v_pts, '.' + color, markersize=8)
    for color, key in [('g', 'green'), ('m', 'magenta')]:
        if points[key]:
            i_pts, v_pts = zip(*points[key])
            plt.plot(i_pts, v_pts, '--' + color, linewidth=4)

    plt.xlabel('$I$')
    plt.ylabel(r'$v_\ast$')
    plt.xlim(-4, 8)
    plt.ylim(minv_c - 5, maxv_c + 5)
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
