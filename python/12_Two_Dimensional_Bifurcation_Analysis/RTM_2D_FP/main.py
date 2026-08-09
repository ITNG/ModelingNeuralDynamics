import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 80.
g_na = 100.
g_l = 0.1
v_k = -100.
v_na = 50.
v_l = -67.


def alpha_m(v):
    return 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))


def alpha_m_p(v):
    num, den = 0.32 * (v + 54), 1 - exp(-(v + 54) / 4)
    num_p, den_p = 0.32, exp(-(v + 54) / 4) / 4
    return (den * num_p - num * den_p) / den ** 2


def beta_m(v):
    return 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)


def beta_m_p(v):
    num, den = 0.28 * (v + 27), exp((v + 27) / 5) - 1
    num_p, den_p = 0.28, exp((v + 27) / 5) / 5
    return (den * num_p - num * den_p) / den ** 2


def alpha_n(v):
    # matlab/12/RTM_2D_FP/alpha_n.m has a removable singularity at v=-52
    # and no guard for it either -- NaN there, same as here, harmless
    # since NaN comparisons are just False in both languages
    with np.errstate(divide='ignore', invalid='ignore'):
        return 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))


def alpha_n_p(v):
    num, den = 0.032 * (v + 52), 1 - exp(-(v + 52) / 5)
    num_p, den_p = 0.032, exp(-(v + 52) / 5) / 5
    return (den * num_p - num * den_p) / den ** 2


def beta_n(v):
    return 0.5 * exp(-(v + 57) / 40)


def beta_n_p(v):
    return -beta_n(v) / 40


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def m_inf_p(v):
    num = (alpha_m(v) + beta_m(v)) * alpha_m_p(v)
    num = num - alpha_m(v) * (alpha_m_p(v) + beta_m_p(v))
    return num / (alpha_m(v) + beta_m(v)) ** 2


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def F(v, i_ext):
    '''dv/dt of the reduced (m=m_inf(v), h=1-n_inf(v)) RTM model'''
    return g_na * m_inf(v) ** 3 * (1 - n_inf(v)) * (v_na - v) \
        + g_k * n_inf(v) ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext


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
        while b_ - a > 1e-10:
            m = (a + b_) / 2
            if F(m, i_ext) * F(a, i_ext) <= 0:
                b_ = m
            else:
                a = m
        roots.append((a + b_) / 2)
    return roots


def jacobian(v, n):
    j00 = g_na * 3 * m_inf(v) ** 2 * m_inf_p(v) * (1 - n) * (v_na - v) \
        - g_na * m_inf(v) ** 3 * (1 - n) - g_k * n ** 4 - g_l
    j01 = -g_na * m_inf(v) ** 3 * (v_na - v) + g_k * 4 * n ** 3 * (v_k - v)
    j10 = alpha_n_p(v) * (1 - n) - beta_n_p(v) * n
    j11 = -alpha_n(v) - beta_n(v)
    return np.array([[j00, j01], [j10, j11]]) / c


i_ext_vec = np.arange(1001) / 1000 * 0.2

points = {'g': [], 'r': [], 'k': [], 'm': [], 'b': []}
for i_ext in i_ext_vec:
    for v in find_fixed_points(i_ext):
        n = n_inf(v)
        e = np.linalg.eigvals(jacobian(v, n))
        if abs(e[0].imag) > 1e-4:
            if e[0].real > 0:
                points['g'].append((i_ext, v))
            if e[0].real < 0:
                points['r'].append((i_ext, v))
        if e[0].real < 0 and e[1].real < 0:
            points['k'].append((i_ext, v))
        if e[0].real * e[1].real < 0:
            points['m'].append((i_ext, v))
        if e[0].real > 0 and e[1].real > 0:
            points['b'].append((i_ext, v))

if __name__ == "__main__":

    plt.figure(figsize=(7, 7))
    for color, pts in points.items():
        if pts:
            i_pts, v_pts = zip(*pts)
            plt.plot(i_pts, v_pts, '.' + color, markersize=3)

    plt.xlim(min(i_ext_vec), max(i_ext_vec))
    plt.ylim(-70, -30)
    plt.xlabel(r'$I$ [$\mu$A/cm$^2$]')
    plt.ylabel(r'$v_\ast$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
