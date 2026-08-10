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

i_ext = 0.

T = 50.  # period of inputs
N = 200
a, b = 0.142, 0.195
g_syn_vec = a + np.arange(N + 1) / N * (b - a)
tau_d = 2.
tau_r = 0.5
tau_peak = tau_r

t_final = 10000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


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


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def tau_peak_function(tau_d, tau_r, tau_d_q):
    dt = 0.01
    dt05 = dt / 2
    s, t = 0., 0.
    s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s * tau_d
    while s_inc > 0:
        t_old, s_inc_old = t, s_inc
        s_tmp = s + dt05 * s_inc
        s_inc_tmp = exp(-(t + dt05) / tau_d_q) * (1 - s_tmp) / tau_r - s_tmp / tau_d
        s = s + dt * s_inc_tmp
        t = t + dt
        s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s / tau_d
    return (t_old * (-s_inc) + t * s_inc_old) / (s_inc_old - s_inc)


def tau_d_q_function(tau_d, tau_r, tau_hat):
    tau_d_q_left = 1.
    while tau_peak_function(tau_d, tau_r, tau_d_q_left) > tau_hat:
        tau_d_q_left /= 2
    tau_d_q_right = tau_r
    while tau_peak_function(tau_d, tau_r, tau_d_q_right) < tau_hat:
        tau_d_q_right *= 2
    while tau_d_q_right - tau_d_q_left > 1e-12:
        tau_d_q_mid = (tau_d_q_left + tau_d_q_right) / 2
        if tau_peak_function(tau_d, tau_r, tau_d_q_mid) <= tau_hat:
            tau_d_q_left = tau_d_q_mid
        else:
            tau_d_q_right = tau_d_q_mid
    return (tau_d_q_left + tau_d_q_right) / 2


tau_dq = tau_d_q_function(tau_d, tau_r, tau_peak)


def spike_times_for(g_syn):
    '''Heun/RK2 integration (plain floats) of a periodically-driven WB
    neuron, returning spike times in the second half of the run (so
    initial transients don't pollute the period statistics)'''
    v = -65.
    m = m_inf(v)
    h = h_inf(v)
    n = n_inf(v)
    q = 0.
    s = 0.
    v_prev = v
    t_spikes = []

    for k in range(m_steps):
        t = k * dt
        if k > 0 and abs(round(t / T) - t / T) < 1e-12:
            q = 1.

        v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v)
                  + g_l * (v_l - v) - g_syn * s * v + i_ext) / c
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
        s_inc = q * (1 - s) / tau_r - s / tau_d
        q_inc = -q / tau_dq

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        s_tmp = s + dt05 * s_inc
        q_tmp = q + dt05 * q_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) - g_syn * s_tmp * v_tmp + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d
        q_inc = -q_tmp / tau_dq

        v_prev = v
        v = v + dt * v_inc
        m = m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc
        s = s + dt * s_inc
        q = q + dt * q_inc

        if v_prev >= -20 and v < -20 and t >= t_final / 2:
            t_next = (k + 1) * dt
            t_spike = (t * (-v - 20) + t_next * (20 + v_prev)) / (v_prev - v)
            t_spikes.append(t_spike)

    return np.array(t_spikes)


n_vec = np.zeros(len(g_syn_vec))
sigma_vec = np.zeros(len(g_syn_vec))
for ijk, g_syn in enumerate(g_syn_vec):
    t_spikes = spike_times_for(g_syn)
    periods = np.diff(t_spikes)
    sigma_vec[ijk] = (periods.max() - periods.min()) / periods.mean()
    n_vec[ijk] = periods.mean() / T

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(g_syn_vec, n_vec, '.k', markersize=5)
    ind = sigma_vec < 1e-3
    plt.plot(g_syn_vec[ind], n_vec[ind], '.r', markersize=10)

    plt.xlim(0.14, 0.19)
    plt.ylim(0, 18)
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\overline{g}_{\rm syn}$')
    plt.ylabel('$n$')
    plt.xticks(np.arange(0.14, 0.191, 0.01))
    plt.yticks(range(5, 16, 5))

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
