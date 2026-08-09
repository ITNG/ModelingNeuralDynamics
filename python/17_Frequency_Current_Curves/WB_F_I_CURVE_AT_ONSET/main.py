import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 9.
g_na = 35.
g_l = 0.1
v_k = -90.
v_na = 55.
v_l = -65.

dt = 0.01
dt05 = dt / 2
N = round(1000 / dt)  # steady-state check window, 1000ms
# matlab's while-loop here has no hard cap (its z=zeros(20000/dt) is just
# an initial buffer size, silently grown if a point needs longer -- some
# points this close to I_c genuinely take >20s of simulated time before
# the 4th spike or the steady-state window criterion fires)
t_max_steps = round(200000 / dt)


def alpha_h(v):
    return 0.35 * exp(-(v + 58) / 20)


def alpha_m(v):
    return 0.1 * (v + 35) / (1 - exp(-(v + 35) / 10))


def alpha_n(v):
    return 0.05 * (v + 34) / (1 - exp(-0.1 * (v + 34)))


def beta_h(v):
    return 5. / (exp(-0.1 * (v + 28)) + 1)


def beta_m(v):
    return 4 * exp(-(v + 60) / 18)


def beta_n(v):
    return 0.625 * exp(-(v + 44) / 80)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def run_to_frequency(i_ext, v, m, h, n):
    '''Heun/RK2 integration (plain floats) of the WB model at fixed i_ext,
    starting from the state left over by the previous i_ext. Returns
    (f, v, m, h, n) -- f=0 if v settles to rest within 20s, else the
    frequency from the 3rd/4th spike interval. 20s (instead of the usual
    5s) because f->0 continuously here (Type I onset, no hysteresis) so
    convergence gets arbitrarily slow right at threshold.'''
    win_maxv = win_minv = v
    win_maxm = win_minm = m
    win_maxh = win_minh = h
    win_maxn = win_minn = n
    num_spikes = 0
    t_spikes = []

    for k in range(1, t_max_steps + 1):
        v_prev = v
        v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext) / c
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc

        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v = v + dt * v_inc
        m = m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc

        win_maxv, win_minv = max(win_maxv, v), min(win_minv, v)
        win_maxm, win_minm = max(win_maxm, m), min(win_minm, m)
        win_maxh, win_minh = max(win_maxh, h), min(win_minh, h)
        win_maxn, win_minn = max(win_maxn, n), min(win_minn, n)

        if v < -20 and v_prev >= -20:
            num_spikes += 1
            t_spikes.append((k * dt * (20 + v_prev) + (k - 1) * dt * (-20 - v))
                             / (v_prev - v))
            if num_spikes == 4:
                return 1000 / (t_spikes[3] - t_spikes[2]), v, m, h, n

        if k % N == 0:
            if (win_maxv - win_minv) < 1e-4 * abs(win_maxv + win_minv) and \
               (win_maxm - win_minm) < 1e-4 * abs(win_maxm + win_minm) and \
               (win_maxh - win_minh) < 1e-4 * abs(win_maxh + win_minh) and \
               (win_maxn - win_minn) < 1e-4 * abs(win_maxn + win_minn):
                return 0., v, m, h, n
            win_maxv = win_minv = v
            win_maxm = win_minm = m
            win_maxh = win_minh = h
            win_maxn = win_minn = n

    raise RuntimeError(f"did not settle within {t_max_steps} steps at I={i_ext}")


i_ext_low, i_ext_high = 0.1600, 0.1602
i_ext_vec = i_ext_low + np.arange(11) / 10 * (i_ext_high - i_ext_low)

f_vec = np.zeros(len(i_ext_vec))
v, m, h, n = -70., m_inf(-70.), 0.7, 0.6
for ijk, i_ext in enumerate(i_ext_vec):
    f_vec[ijk], v, m, h, n = run_to_frequency(i_ext, v, m, h, n)

# fit f = C*sqrt(I-I_c) near onset: search over candidate I_c between the
# last non-firing and first firing i_ext for the flattest C = f/sqrt(I-I_c)
ind = np.where(f_vec > 0)[0]
i_ext_vec_0 = i_ext_vec[ind]
f_vec_0 = f_vec[ind]
I_c_low = i_ext_vec[ind.min() - 1]
I_c_high = i_ext_vec[ind.min()]

alpha_vec = np.arange(101) / 100
C_vec = np.zeros(len(alpha_vec))
err_vec = np.zeros(len(alpha_vec))
for ijk, alpha in enumerate(alpha_vec):
    I_c_trial = I_c_low * alpha + I_c_high * (1 - alpha)
    # alpha=0 sets I_c_trial exactly to the first spiking i_ext -> 1/0
    # there; matlab's min() silently skips that NaN entry, so match that
    with np.errstate(divide='ignore', invalid='ignore'):
        y = f_vec_0 / np.sqrt(i_ext_vec_0 - I_c_trial)
        C_vec[ijk] = y.mean()
        err_vec[ijk] = (y.max() - y.min()) / y.mean()

ind_best = np.where(err_vec == np.nanmin(err_vec))[0].min()
alpha = alpha_vec[ind_best]
I_c = I_c_low * alpha + I_c_high * (1 - alpha)
C = C_vec[ind_best]

if __name__ == "__main__":

    print(f"I_c = {I_c}")
    print(f"C = {C}")

    plt.figure(figsize=(7, 4))
    I = I_c + np.arange(1001) / 1000 * (i_ext_high - I_c)
    plt.plot(I, C * np.sqrt(I - I_c), '-r', linewidth=2)
    plt.plot([i_ext_low, I_c], [0, 0], '-r', linewidth=4)
    plt.plot(i_ext_vec, f_vec, '.k', markersize=20)
    plt.xlabel('$I$')
    plt.ylabel('$f$')
    plt.xlim(i_ext_low, i_ext_high)
    plt.ylim(0, f_vec.max() * 1.1)
    plt.xticks([0.1600, 0.1601, 0.1602])
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
