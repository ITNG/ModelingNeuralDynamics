import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 224.
g_na = 112.
g_l = 0.5
v_k = -90.
v_na = 60.
v_l = -70.


def alpha_m(v):
    return 40 * (75.5 - v) / (exp((75.5 - v) / 13.5) - 1)


def alpha_m_p(v):
    num, den = 40 * (75.5 - v), exp((75.5 - v) / 13.5) - 1
    num_p, den_p = -40., -exp((75.5 - v) / 13.5) / 13.5
    return (den * num_p - num * den_p) / den ** 2


def beta_m(v):
    return 1.2262 / exp(v / 42.248)


def beta_m_p(v):
    return -beta_m(v) / 42.248


def alpha_n(v):
    return (95 - v) / (exp((95 - v) / 11.8) - 1)


def alpha_n_p(v):
    num, den = 95 - v, exp((95 - v) / 11.8) - 1
    num_p, den_p = -1., -(den + 1) / 11.8
    return (den * num_p - num * den_p) / den ** 2


def beta_n(v):
    return 0.025 / exp(v / 22.222)


def beta_n_p(v):
    return -beta_n(v) / 22.222


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def m_inf_p(v):
    num = (alpha_m(v) + beta_m(v)) * alpha_m_p(v)
    num = num - alpha_m(v) * (alpha_m_p(v) + beta_m_p(v))
    return num / (alpha_m(v) + beta_m(v)) ** 2


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def F(v, i_ext):
    '''dv/dt of the reduced (m=m_inf(v), h=0.36-n_inf(v)) Erisir model'''
    return g_na * m_inf(v) ** 3 * (0.36 - n_inf(v)) * (v_na - v) \
        + g_k * n_inf(v) ** 2 * (v_k - v) + g_l * (v_l - v) + i_ext


def find_fixed_points(i_ext):
    '''bracket sign changes of F on a v-grid, then bisect each one'''
    v_min = min(v_k, v_l + i_ext / g_l)
    v_max = max(v_na, v_l + i_ext / g_l)
    n_v = 1000
    v_vec = v_min + np.arange(n_v + 1) / n_v * (v_max - v_min)
    f_vec = F(v_vec, i_ext)

    roots = []
    for i in np.where(f_vec[:-1] * f_vec[1:] <= 0)[0]:
        a, b_ = v_vec[i], v_vec[i + 1]
        while b_ - a > 1e-12:
            m = (a + b_) / 2
            if F(m, i_ext) * F(a, i_ext) <= 0:
                b_ = m
            else:
                a = m
        roots.append((a + b_) / 2)
    return roots


def jacobian(v, n):
    j00 = g_na * 3 * m_inf(v) ** 2 * m_inf_p(v) * (0.36 - n) * (v_na - v) \
        - g_na * m_inf(v) ** 3 * (0.36 - n) - g_k * n ** 2 - g_l
    j01 = -g_na * m_inf(v) ** 3 * (v_na - v) + g_k * 2 * n * (v_k - v)
    j10 = alpha_n_p(v) * (1 - n) - beta_n_p(v) * n
    j11 = -alpha_n(v) - beta_n(v)
    return np.array([[j00, j01], [j10, j11]]) / c


i_ext_vec = np.arange(1001) / 1000 * 7

points = {'g': [], 'r': [], 'k': [], 'm': [], 'b': []}
for i_ext in i_ext_vec:
    for v in find_fixed_points(i_ext):
        n = n_inf(v)
        e = np.linalg.eigvals(jacobian(v, n))
        if abs(e[0].imag) > 1e-6:
            if e[0].real > 0:
                points['g'].append((i_ext, v))
            if e[0].real < 0:
                points['r'].append((i_ext, v))
        if abs(e[0].imag) < 1e-6:
            if e[0].real < 0 and e[1].real < 0:
                points['k'].append((i_ext, v))
            if e[0].real * e[1].real <= 0:
                points['m'].append((i_ext, v))
            if e[0].real > 0 and e[1].real > 0:
                points['b'].append((i_ext, v))

if __name__ == "__main__":

    plt.figure(figsize=(7, 7))
    for color, pts in points.items():
        if not pts:
            continue
        i_pts, v_pts = zip(*pts)
        if color in ('m', 'b'):
            # matlab plots these two as connected dashed lines, not dots
            plt.plot(i_pts, v_pts, '--' + color, linewidth=2)
        else:
            plt.plot(i_pts, v_pts, '.' + color, markersize=3)

    plt.xlim(min(i_ext_vec), max(i_ext_vec))
    plt.ylim(-90, 0)
    plt.xlabel(r'$I$ [$\mu$A/cm$^2$]')
    plt.ylabel(r'$v_\ast$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
