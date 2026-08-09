import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

c = 1.
g_k = 36.
g_na = 120.
g_l = 0.3
v_k = -82.
v_na = 45.
v_l = -59.


def alpha_m(v):
    with np.errstate(divide='ignore', invalid='ignore'):
        out = (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))
    return out if abs(v + 45) > 1e-8 else 1.0


def alpha_m_p(v):
    num, den = (v + 45) / 10, 1 - exp(-(v + 45) / 10)
    num_p, den_p = 1 / 10, exp(-(v + 45) / 10) / 10
    return (den * num_p - num * den_p) / den ** 2


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_m_p(v):
    return -(4 / 18) * exp(-(v + 70) / 18)


def alpha_n(v):
    return 0.01 * (-60.0 - v) / (exp((-60 - v) / 10) - 1)


def alpha_n_p(v):
    num, den = 0.01 * (-60.0 - v), exp((-60 - v) / 10) - 1
    num_p, den_p = -0.01, -(den + 1) * 0.1
    return (den * num_p - num * den_p) / den ** 2


def beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def beta_n_p(v):
    return -beta_n(v) / 80


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def m_inf_p(v):
    num = (alpha_m(v) + beta_m(v)) * alpha_m_p(v)
    num = num - alpha_m(v) * (alpha_m_p(v) + beta_m_p(v))
    return num / (alpha_m(v) + beta_m(v)) ** 2


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def f(v):
    '''dv/dt at i_ext=0 of the reduced (m=m_inf(v), h=0.83-n_inf(v)) HH model'''
    return g_na * m_inf(v) ** 3 * (0.83 - n_inf(v)) * (v_na - v) \
        + g_k * n_inf(v) ** 4 * (v_k - v) + g_l * (v_l - v)


def find_fixed_point(i_ext):
    v_left, v_right = -100., 50.
    while v_right - v_left > 1e-10:
        v_c = (v_left + v_right) / 2
        if (f(v_c) + i_ext) * (f(v_left) + i_ext) > 0:
            v_left = v_c
        else:
            v_right = v_c
    return (v_left + v_right) / 2


def jacobian(v_c, n_c):
    j00 = g_na * 3 * m_inf(v_c) ** 2 * m_inf_p(v_c) * (0.83 - n_c) * (v_na - v_c) \
        - g_na * m_inf(v_c) ** 3 * (0.83 - n_c) - g_k * n_c ** 4 - g_l
    j01 = -g_na * m_inf(v_c) ** 3 * (v_na - v_c) + 4 * g_k * n_c ** 3 * (v_k - v_c)
    j10 = alpha_n_p(v_c) * (1 - n_c) - beta_n_p(v_c) * n_c
    j11 = -alpha_n(v_c) - beta_n(v_c)
    return np.array([[j00, j01], [j10, j11]]) / c


i_ext_vec = np.arange(1001) / 1000 * 15

real_part = np.zeros(len(i_ext_vec))
imag_part = np.zeros(len(i_ext_vec))
for ijk, i_ext in enumerate(i_ext_vec):
    v_c = find_fixed_point(i_ext)
    n_c = n_inf(v_c)
    e = np.linalg.eigvals(jacobian(v_c, n_c))
    real_part[ijk] = e[0].real
    imag_part[ijk] = e[0].imag

if __name__ == "__main__":

    A, B, C, D = -1, 1, -1, 1
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(real_part, imag_part, color='k', linewidth=2)
    ax.plot(real_part, -imag_part, color='k', linewidth=2)

    ind = np.where((real_part[:-1] < 0.1) & (real_part[1:] > 0.1))[0][0]
    x, y = real_part[ind], imag_part[ind]
    v = np.array([real_part[ind] - real_part[ind - 200], imag_part[ind] - imag_part[ind - 200]])
    draw_arrow(ax, [A, B], [C, D], x, y, v, epsilon=0.1, width=2)
    draw_arrow(ax, [A, B], [C, D], x, -y, [v[0], -v[1]], epsilon=0.1, width=2)

    ax.plot([0, 0], [-1, 1], color='b', linewidth=1)
    ax.set_xlim(A, B)
    ax.set_ylim(C, D)
    ax.set_box_aspect(1)
    ax.set_xlabel(r'${\rm Re}\,(\lambda)$')
    ax.set_ylabel(r'${\rm Im}\,(\lambda)$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
