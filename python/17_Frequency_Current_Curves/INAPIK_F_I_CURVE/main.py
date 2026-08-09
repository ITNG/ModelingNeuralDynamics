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

dt = 0.002
dt05 = dt / 2
N = round(1000 / dt)  # steady-state check window, 1000ms
# matlab's while-loop has no hard cap (z=zeros(5000/dt) is just an
# initial buffer, silently grown further if a point needs longer)
t_max_steps = round(50000 / dt)


def m_inf(v):
    return 1. / (1 + exp((-20 - v) / 15))


def n_inf(v):
    return 1. / (1 + exp((-25 - v) / 5))


def run_to_frequency(i_ext, v, m, n):
    '''Heun/RK2 integration (plain floats) of the INaP+IK model at fixed
    i_ext, starting from the state left over by the previous i_ext.
    Returns (f, v, m, n) -- f=0 if v settles to rest, else the frequency
    from the 3rd/4th spike interval.'''
    win_maxv = win_minv = v
    win_maxm = win_minm = m
    win_maxn = win_minn = n
    num_spikes = 0
    t_spikes = []

    for k in range(1, t_max_steps + 1):
        v_prev = v
        v_inc = (g_na * m * (v_na - v) + g_k * n * (v_k - v) + g_l * (v_l - v) + i_ext) / c
        n_inc = (n_inf(v) - n) / tau_n

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        n_tmp = n + dt05 * n_inc

        v_inc = (g_na * m_tmp * (v_na - v_tmp) + g_k * n_tmp * (v_k - v_tmp) + g_l * (v_l - v_tmp) + i_ext) / c
        n_inc = (n_inf(v_tmp) - n_tmp) / tau_n

        v = v + dt * v_inc
        m = m_inf(v)
        n = n + dt * n_inc

        win_maxv, win_minv = max(win_maxv, v), min(win_minv, v)
        win_maxm, win_minm = max(win_maxm, m), min(win_minm, m)
        win_maxn, win_minn = max(win_maxn, n), min(win_minn, n)

        if v < -20 and v_prev >= -20:
            num_spikes += 1
            t_spikes.append((k * dt * (20 + v_prev) + (k - 1) * dt * (-20 - v))
                             / (v_prev - v))
            if num_spikes == 4:
                return 1000 / (t_spikes[3] - t_spikes[2]), v, m, n

        if k % N == 0:
            if (win_maxv - win_minv) < 1e-4 * abs(win_maxv + win_minv) and \
               (win_maxm - win_minm) < 1e-4 * abs(win_maxm + win_minm) and \
               (win_maxn - win_minn) < 1e-4 * abs(win_maxn + win_minn):
                return 0., v, m, n
            win_maxv = win_minv = v
            win_maxm = win_minm = m
            win_maxn = win_minn = n

    raise RuntimeError(f"did not settle within {t_max_steps} steps at I={i_ext}")


i_ext_vec = -4 + np.arange(31) / 30 * 12

f_forward = np.zeros(len(i_ext_vec))
v, m, n = -70., m_inf(-70.), 0.6
for ijk, i_ext in enumerate(i_ext_vec):
    f_forward[ijk], v, m, n = run_to_frequency(i_ext, v, m, n)

f_backward = np.zeros(len(i_ext_vec))
for ijk in range(len(i_ext_vec) - 1, -1, -1):
    f_backward[ijk], v, m, n = run_to_frequency(i_ext_vec[ijk], v, m, n)

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
    plt.xlim(-3, 7)
    plt.ylim(0, 1200)
    plt.xlabel('$I$')
    plt.ylabel('$f$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
