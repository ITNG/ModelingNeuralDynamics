import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

alpha = 1.  # sharpness of the input pulses
Period = 25.  # period of the input pulses
N = 2000
ave = np.mean(np.exp(alpha * np.cos(np.pi * np.arange(N) / N) ** 2) - 1)


def shape(t):
    '''the input pulses always have the same shape, with temporal
    average 1; amplitude is applied by the caller.'''
    return (np.exp(alpha * np.cos(np.pi * t / Period) ** 2) - 1) / ave


c = 1.
g_k, g_na, g_l = 80., 100., 0.2
v_k, v_na, v_l = -100., 50., -67.

i_ext_low, i_ext_high = 0., 2.
i_ext_vec = i_ext_low + np.arange(201) / 200 * (i_ext_high - i_ext_low)

dt = 0.01
dt05 = dt / 2
t_final = 1000.
m_steps = round(t_final / dt)
t = np.arange(m_steps + 1) * dt
shape_store = shape(t)


def m_inf(v):
    alpha_m = 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))
    beta_m = 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)
    return alpha_m / (alpha_m + beta_m)


def alpha_h(v):
    return 0.128 * exp(-(v + 50) / 18)


def beta_h(v):
    return 4. / (1 + exp(-(v + 27) / 5))


def alpha_n(v):
    return 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))


def beta_n(v):
    return 0.5 * exp(-(v + 57) / 40)


def step(v, m, h, n, i_ext_old, i_ext_mid):
    v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext_old) / c
    h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
    n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n

    v_tmp = v + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h + dt05 * h_inc
    n_tmp = n + dt05 * n_inc

    v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
              + g_l * (v_l - v_tmp) + i_ext_mid) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

    v_new = v + dt * v_inc
    m_new = m_inf(v_new)
    h_new = h + dt * h_inc
    n_new = n + dt * n_inc
    return v_new, m_new, h_new, n_new


def f_i_curve_constant():
    '''Frequency from the interval between the 1st and 2nd spikes under
    a constant drive, continuing from the previous drive's final state;
    if two spikes don't occur within t_final, frequency is taken as 0.'''
    v, m, h, n = -70., m_inf(-70.), 0.7, 0.6
    f_vec = np.zeros(len(i_ext_vec))

    for ijk, i_ext in enumerate(i_ext_vec):
        t_spikes = []
        v_old = v
        f = 0.
        for k in range(1, m_steps + 1):
            v, m, h, n = step(v, m, h, n, i_ext, i_ext)
            if v < -20 and v_old >= -20:
                t_spikes.append((k * dt * (20 + v_old) + (k - 1) * dt * (-20 - v)) / (v_old - v))
                if len(t_spikes) == 2:
                    f = 1000. / (t_spikes[1] - t_spikes[0])
                    break
            v_old = v
        f_vec[ijk] = f

    return f_vec


def f_i_curve_pulsed():
    '''Firing rate (spike count over the fixed t_final window) under a
    periodic pulsed drive of the given time-average amplitude,
    continuing from the previous drive's final state.'''
    v, m, h, n = -70., m_inf(-70.), 0.7, 0.6
    f_vec = np.zeros(len(i_ext_vec))

    for ijk, i_ext in enumerate(i_ext_vec):
        i_ext_p = i_ext * shape_store
        num_spikes = 0
        v_old = v
        for k in range(1, m_steps + 1):
            v, m, h, n = step(v, m, h, n, i_ext_p[k - 1], (i_ext_p[k - 1] + i_ext_p[k]) / 2)
            if v < -20 and v_old >= -20:
                num_spikes += 1
            v_old = v
        f_vec[ijk] = num_spikes / t_final * 1000

    return f_vec


f_vec_constant = f_i_curve_constant()
f_vec_pulsed = f_i_curve_pulsed()

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(i_ext_vec, f_vec_constant, '.r', markersize=8, label='constant drive')
    ax.plot(i_ext_vec, f_vec_pulsed, '.b', markersize=8, label='pulsed drive')
    ax.set_xlabel(r'$I$ [$\mu$A/cm$^2$]')
    ax.set_ylabel('$f$')
    ax.legend()

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
