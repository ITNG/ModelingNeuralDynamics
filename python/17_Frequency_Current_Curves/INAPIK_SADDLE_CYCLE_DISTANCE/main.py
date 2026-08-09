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


w_grid = -100 + np.arange(10001) / 10000 * 150
f_grid = f(w_grid)


def find_saddle(I):
    '''bracket sign changes of f+I on w_grid, bisect each, and return the
    (v, n) of whichever fixed point is a saddle'''
    fi = f_grid + I
    v_star = n_star = None
    for j in np.where(fi[:-1] * fi[1:] <= 0)[0]:
        w_low, w_high = w_grid[j], w_grid[j + 1]
        while w_high - w_low > 1e-12:
            w_c = (w_low + w_high) / 2
            if (f(w_c) + I) * (f(w_high) + I) <= 0:
                w_low = w_c
            else:
                w_high = w_c
        v_c = (w_low + w_high) / 2
        n_c = n_inf(v_c)
        j00 = g_na * m_inf_p(v_c) * (v_na - v_c) - g_na * m_inf(v_c) - g_k * n_c - g_l
        j01 = g_k * (v_k - v_c)
        j10 = n_inf_p(v_c) / tau_n
        j11 = -1 / tau_n
        e = np.linalg.eigvals(np.array([[j00, j01], [j10, j11]]))
        if abs(e[0].imag) < 1e-12 and e[0].real * e[1].real < 0:
            v_star, n_star = v_c, n_c
    return v_star, n_star


t_final = 20.
dt = 0.001
dt05 = dt / 2
m_steps = round(t_final / dt)

low, high = -1.38, 0.
i_ext_vec = low + np.arange(101) / 100 * (high - low)

d_vec = np.zeros(len(i_ext_vec))
for ijk, i_ext in enumerate(i_ext_vec):
    v = np.zeros(m_steps + 1)
    n = np.zeros(m_steps + 1)
    v[0], n[0] = -30., 0.2

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

    v_star, n_star = find_saddle(i_ext)
    d_vec[ijk] = np.sqrt(((v - v_star) ** 2 + (n - n_star) ** 2).min())

if __name__ == "__main__":

    plt.figure(figsize=(7, 3.5))
    plt.plot(i_ext_vec, d_vec, '-k', linewidth=2)
    plt.xlabel('$I$')
    plt.ylabel('$d$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
