import numpy as np
from numpy import exp

c = 1.
g_k = 36.
g_na = 120.
g_l = 0.3
v_k = -82.
v_na = 45.
v_l = -59.


def alpha_m(v):
    # removable singularity at v=-45 (0/0) -> L'Hopital limit is 1
    with np.errstate(divide='ignore', invalid='ignore'):
        out = (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))
    return np.where(np.abs(v + 45) > 1e-8, out, 1.0)


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


v = -100 + np.arange(10001) / 10000 * 150

i_ext_vec = np.arange(1001) / 1000 * 15
num_fp = np.zeros(len(i_ext_vec), dtype=int)

for ijk, i_ext in enumerate(i_ext_vec):
    f = g_na * m_inf(v) ** 3 * (0.83 - n_inf(v)) * (v_na - v) \
        + g_k * n_inf(v) ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext
    num_fp[ijk] = np.sum(f[:-1] * f[1:] <= 0)

if __name__ == "__main__":
    print(num_fp.max())
    print(num_fp.min())
