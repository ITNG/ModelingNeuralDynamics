import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

# the backward (unstable-cycle) integration below genuinely blows up for a
# handful of i_ext values before the "if Bmax<200" check discards it --
# same as matlab, which silently lets Inf/NaN propagate there too
np.seterr(all='ignore')

c = 1.
g_k = 36.
g_na = 120.
g_l = 0.3
v_k = -82.
v_na = 45.
v_l = -59.


def alpha_h(v):
    return 0.07 * exp(-(v + 70) / 20)


def alpha_m(v):
    with np.errstate(divide='ignore', invalid='ignore'):
        out = (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))
    return out if abs(v + 45) > 1e-8 else 1.0


def alpha_m_p(v):
    num, den = (v + 45) / 10, 1 - exp(-(v + 45) / 10)
    num_p, den_p = 1 / 10, exp(-(v + 45) / 10) / 10
    return (den * num_p - num * den_p) / den ** 2


def alpha_n(v):
    return 0.01 * (-60. - v) / (exp((-60 - v) / 10) - 1)


def alpha_n_p(v):
    num, den = 0.01 * (-60. - v), exp((-60 - v) / 10) - 1
    num_p, den_p = -0.01, -(den + 1) * 0.1
    return (den * num_p - num * den_p) / den ** 2


def beta_h(v):
    return 1. / (exp(-(v + 40) / 10) + 1)


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_m_p(v):
    return -(4 / 18) * exp(-(v + 70) / 18)


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
    '''dv/dt of the reduced (m=m_inf(v), h=0.83-n_inf(v)) HH model, i_ext=0'''
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


def jacobian(v, n):
    j00 = g_na * 3 * m_inf(v) ** 2 * m_inf_p(v) * (0.83 - n) * (v_na - v) \
        - g_na * m_inf(v) ** 3 * (0.83 - n) - g_k * n ** 4 - g_l
    j01 = -g_na * m_inf(v) ** 3 * (v_na - v) + 4 * g_k * n ** 3 * (v_k - v)
    j10 = alpha_n_p(v) * (1 - n) - beta_n_p(v) * n
    j11 = -alpha_n(v) - beta_n(v)
    return np.array([[j00, j01], [j10, j11]]) / c


t_final = 300.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)
tail_start = round(2 / 3 * m_steps)


def envelope(v0, n0, i_ext, direction):
    '''Heun integration (plain floats), returns max/min of v over the
    tail 1/3 of the run. direction=+1 forward (settle onto the stable
    cycle), direction=-1 backward (trace the unstable cycle out from
    just next to the stable fixed point).'''
    v, n = v0, n0
    vmax, vmin = -np.inf, np.inf
    for k in range(m_steps):
        if not np.isfinite(v):
            # diverged (this happens for a few i_ext when integrating the
            # unstable cycle backward) -- report as blown up, same as
            # matlab's Inf/NaN propagating through to fail the Bmax<200 check
            return np.inf, -np.inf
        m = m_inf(v)
        h = 0.83 - n
        v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v)
                  + g_l * (v_l - v) + i_ext) / c
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
        v_tmp = v + direction * dt05 * v_inc
        n_tmp = n + direction * dt05 * n_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = 0.83 - n_tmp
        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        v = v + direction * dt * v_inc
        n = n + direction * dt * n_inc
        if k >= tail_start:
            if v > vmax:
                vmax = v
            if v < vmin:
                vmin = v
    return vmax, vmin


i_ext_vec = 5 + np.arange(101) / 100 * (10 - 5)

green = []  # unstable fixed point (unstable spiral)
red = []    # stable fixed point (stable spiral)
uc = []     # (i_ext, Bmax, Bmin) unstable-cycle envelope, from near the red branch
sc = []     # (i_ext, Amax, Amin) stable-cycle envelope

for i_ext in i_ext_vec:
    v_c = find_fixed_point(i_ext)
    n_c = n_inf(v_c)
    e = np.linalg.eigvals(jacobian(v_c, n_c))

    if e[0].real < 0 and e[1].real < 0 and abs(e[0].imag) > 1e-4:
        red.append((i_ext, v_c))
        # start near the stable fixed point, but not quite, and run
        # backward in time to find the unstable cycle
        Bmax, Bmin = envelope(v_c * 1.01, n_c * 0.99, i_ext, direction=-1)
        if Bmax < 200:
            uc.append((i_ext, Bmax, Bmin))

    if e[0].real > 0 and e[1].real > 0 and abs(e[0].imag) > 1e-4:
        green.append((i_ext, v_c))

    # start far away to settle onto the stable cycle
    Amax, Amin = envelope(70., 0.5, i_ext, direction=1)
    if Amax - Amin > 2:
        sc.append((i_ext, Amax, Amin))

if __name__ == "__main__":

    plt.figure(figsize=(7, 7))

    if green:
        i_pts, v_pts = zip(*green)
        plt.plot(i_pts, v_pts, '--g', linewidth=4)
    if red:
        i_pts, v_pts = zip(*red)
        plt.plot(i_pts, v_pts, '-r', linewidth=4)

    i_sc, amax, amin = zip(*sc)
    plt.plot(i_sc, amax, '-k', linewidth=1)
    plt.plot(i_sc, amin, '-k', linewidth=1)

    i_uc, bmax, bmin = zip(*uc)
    i_uc = [i_sc[0], *i_uc]
    bmax = [amax[0], *bmax]
    bmin = [amin[0], *bmin]
    plt.plot(i_uc, bmax, '--k', linewidth=2)
    plt.plot(i_uc, bmin, '--k', linewidth=2)

    plt.xticks([])
    plt.ylabel('$v$')
    plt.xlim(5, 10)
    plt.ylim(-90, 50)
    plt.text(4.5, -100, r'$I_\ast \approx 5.25$', fontsize=14)
    plt.text(7.4 - 0.75, -100, r'$I_c \approx 7.4$', fontsize=14)
    plt.plot([5.25, 5.25], [-92, -87], '-k', linewidth=1)
    plt.plot([7.4, 7.4], [-92, -87], '-k', linewidth=1)
    plt.gca().set_box_aspect(1)
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
