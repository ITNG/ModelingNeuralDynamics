import itertools

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

g_h = 1.
v_h = -32.9

dt = 0.01
dt05 = dt / 2
N = round(1000 / dt)  # steady-state check window, 1000ms


def alpha_h(v):
    return 0.128 * exp(-(v + 50) / 18)


def alpha_m(v):
    return 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))


def alpha_n(v):
    return 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))


def beta_h(v):
    return 4. / (1 + exp(-(v + 27) / 5))


def beta_m(v):
    return 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)


def beta_n(v):
    return 0.5 * exp(-(v + 57) / 40)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def r_inf(v):
    return 1. / (1 + exp((v + 84) / 10.2))


def tau_r(v):
    return 1. / (exp(-14.59 - 0.086 * v) + exp(-1.87 + 0.0701 * v))


def run_to_frequency(i_ext, v, m, h, n, r):
    '''Heun/RK2 integration (plain floats) of the RTM model with an added
    I_h current, at fixed i_ext, starting from the state left over by the
    previous i_ext. Returns (f, v, m, h, n, r) -- f=0 if v settles to rest,
    else the frequency from the 3rd/4th spike interval.'''
    win_maxv = win_minv = v
    win_maxm = win_minm = m
    win_maxh = win_minh = h
    win_maxn = win_minn = n
    win_maxr = win_minr = r
    num_spikes = 0
    t_spikes = []

    for k in itertools.count(1):
        v_prev = v
        v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v)
                  + g_h * r * (v_h - v) + i_ext) / c
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
        r_inc = (r_inf(v) - r) / tau_r(v)

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        r_tmp = r + dt05 * r_inc

        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_h * r_tmp * (v_h - v_tmp) + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        r_inc = (r_inf(v_tmp) - r_tmp) / tau_r(v_tmp)

        v = v + dt * v_inc
        m = m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc
        r = r + dt * r_inc

        win_maxv, win_minv = max(win_maxv, v), min(win_minv, v)
        win_maxm, win_minm = max(win_maxm, m), min(win_minm, m)
        win_maxh, win_minh = max(win_maxh, h), min(win_minh, h)
        win_maxn, win_minn = max(win_maxn, n), min(win_minn, n)
        win_maxr, win_minr = max(win_maxr, r), min(win_minr, r)

        if v < -20 and v_prev >= -20:
            num_spikes += 1
            t_spikes.append((k * dt * (20 + v_prev) + (k - 1) * dt * (-20 - v))
                             / (v_prev - v))
            if num_spikes == 4:
                return 1000 / (t_spikes[3] - t_spikes[2]), v, m, h, n, r

        if k % N == 0:
            if (win_maxv - win_minv) < 1e-4 * abs(win_maxv + win_minv) and \
               (win_maxm - win_minm) < 1e-4 * abs(win_maxm + win_minm) and \
               (win_maxh - win_minh) < 1e-4 * abs(win_maxh + win_minh) and \
               (win_maxn - win_minn) < 1e-4 * abs(win_maxn + win_minn) and \
               (win_maxr - win_minr) < 1e-4 * abs(win_maxr + win_minr):
                return 0., v, m, h, n, r
            win_maxv = win_minv = v
            win_maxm = win_minm = m
            win_maxh = win_minh = h
            win_maxn = win_minn = n
            win_maxr = win_minr = r


i_ext_low, i_ext_high = -3.2, -3.185
i_ext_vec = i_ext_low + np.arange(31) / 30 * (i_ext_high - i_ext_low)

f_forward = np.zeros(len(i_ext_vec))
v, m, h, n, r = -70., m_inf(-70.), 0.7, 0.6, 0.2
for ijk, i_ext in enumerate(i_ext_vec):
    f_forward[ijk], v, m, h, n, r = run_to_frequency(i_ext, v, m, h, n, r)

f_backward = np.zeros(len(i_ext_vec))
for ijk in range(len(i_ext_vec) - 1, -1, -1):
    f_backward[ijk], v, m, h, n, r = run_to_frequency(i_ext_vec[ijk], v, m, h, n, r)

ind = np.where(f_forward == 0)[0].max()
I_c = (i_ext_vec[ind] + i_ext_vec[ind + 1]) / 2

ind = np.where(f_backward == 0)[0].max()
I_star = (i_ext_vec[ind] + i_ext_vec[ind + 1]) / 2

if __name__ == "__main__":

    print(f"I_c = {I_c}")
    print(f"I_star = {I_star}")

    plt.figure(figsize=(7, 3.5))
    plt.plot(i_ext_vec, f_forward, '.k', markersize=15, label='forward')
    plt.plot(i_ext_vec, f_backward, 'ok', markersize=10, markerfacecolor='none',
              linewidth=1, label='backward')
    plt.xlim(i_ext_low, i_ext_high)
    plt.ylim(0, max(f_forward.max(), f_backward.max()) * 1.1)
    plt.xlabel('$I$')
    plt.ylabel('$f$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
