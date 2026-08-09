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

dt = 0.01
dt05 = dt / 2
N = round(1000 / dt)  # steady-state check window, 1000ms
t_max_steps = round(5000 / dt)


def alpha_m(v):
    with np.errstate(divide='ignore', invalid='ignore'):
        out = (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))
    return out if abs(v + 45) > 1e-8 else 1.0


def alpha_n(v):
    return 0.01 * (-60. - v) / (exp((-60 - v) / 10) - 1)


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def run_to_frequency(i_ext, v, n):
    '''Heun/RK2 integration (plain floats) of the reduced HH model
    (m=m_inf(v), h=0.83-n) at fixed i_ext, starting from the state left
    over by the previous i_ext. Returns (f, v, n) -- f=0 if v settles to
    rest within 5s, else the frequency from the 3rd/4th spike interval.'''
    m = m_inf(v)
    h = 0.83 - n
    win_maxv = win_minv = v
    win_maxm = win_minm = m
    win_maxh = win_minh = h
    win_maxn = win_minn = n
    num_spikes = 0
    t_spikes = []

    for k in range(1, t_max_steps + 1):
        v_prev = v
        v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext) / c
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n

        v_tmp = v + dt05 * v_inc
        n_tmp = n + dt05 * n_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = 0.83 - n_tmp

        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v = v + dt * v_inc
        n = n + dt * n_inc
        m = m_inf(v)
        h = 0.83 - n

        win_maxv, win_minv = max(win_maxv, v), min(win_minv, v)
        win_maxm, win_minm = max(win_maxm, m), min(win_minm, m)
        win_maxh, win_minh = max(win_maxh, h), min(win_minh, h)
        win_maxn, win_minn = max(win_maxn, n), min(win_minn, n)

        if v < -20 and v_prev >= -20:
            num_spikes += 1
            t_spikes.append((k * dt * (20 + v_prev) + (k - 1) * dt * (-20 - v))
                             / (v_prev - v))
            if num_spikes == 4:
                return 1000 / (t_spikes[3] - t_spikes[2]), v, n

        if k % N == 0:
            if (win_maxv - win_minv) < 1e-4 * abs(win_maxv + win_minv) and \
               (win_maxm - win_minm) < 1e-4 * abs(win_maxm + win_minm) and \
               (win_maxh - win_minh) < 1e-4 * abs(win_maxh + win_minh) and \
               (win_maxn - win_minn) < 1e-4 * abs(win_maxn + win_minn):
                return 0., v, n
            win_maxv = win_minv = v
            win_maxm = win_minm = m
            win_maxh = win_minh = h
            win_maxn = win_minn = n

    raise RuntimeError(f"did not settle within {t_max_steps} steps at I={i_ext}")


i_ext_vec = 3 + np.arange(31) / 30 * 10

f_forward = np.zeros(len(i_ext_vec))
v, n = -70., 0.6
for ijk, i_ext in enumerate(i_ext_vec):
    f_forward[ijk], v, n = run_to_frequency(i_ext, v, n)

f_backward = np.zeros(len(i_ext_vec))
for ijk in range(len(i_ext_vec) - 1, -1, -1):
    f_backward[ijk], v, n = run_to_frequency(i_ext_vec[ijk], v, n)

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
    plt.xlim(3, 13)
    plt.ylim(0, 100)
    plt.xlabel('$I$')
    plt.ylabel('$f$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
