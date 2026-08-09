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

T = 50.  # input period
tau_d = 2.
tau_r = 0.5
tau_peak = tau_r

dt = 0.01
dt05 = dt / 2


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


def simulate(g_syn, t_final):
    '''Heun/RK2 integration (numpy arrays) of a WB neuron periodically
    driven by an excitatory synapse (q reset to 1 every T ms)'''
    m_steps = round(t_final / dt)
    v = np.zeros(m_steps + 1)
    m = np.zeros(m_steps + 1)
    h = np.zeros(m_steps + 1)
    n = np.zeros(m_steps + 1)
    q = np.zeros(m_steps + 1)
    s = np.zeros(m_steps + 1)
    v[0] = -65.
    m[0] = m_inf(v[0])
    h[0] = h_inf(v[0])
    n[0] = n_inf(v[0])

    for k in range(m_steps):
        t = k * dt
        if k > 0 and abs(round(t / T) - t / T) < 1e-12:
            q[k] = 1.

        v_inc = (g_k * n[k] ** 4 * (v_k - v[k]) + g_na * m[k] ** 3 * h[k] * (v_na - v[k])
                  + g_l * (v_l - v[k]) - g_syn * s[k] * v[k] + i_ext) / c
        h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]
        n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
        s_inc = q[k] * (1 - s[k]) / tau_r - s[k] / tau_d
        q_inc = -q[k] / tau_dq

        v_tmp = v[k] + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h[k] + dt05 * h_inc
        n_tmp = n[k] + dt05 * n_inc
        s_tmp = s[k] + dt05 * s_inc
        q_tmp = q[k] + dt05 * q_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) - g_syn * s_tmp * v_tmp + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d
        q_inc = -q_tmp / tau_dq

        v[k + 1] = v[k] + dt * v_inc
        m[k + 1] = m_inf(v[k + 1])
        h[k + 1] = h[k] + dt * h_inc
        n[k + 1] = n[k] + dt * n_inc
        s[k + 1] = s[k] + dt * s_inc
        q[k + 1] = q[k] + dt * q_inc

    return v


def spike_times(v, t_final):
    m_steps = round(t_final / dt)
    t = np.arange(m_steps + 1) * dt
    ind = np.where((v[:-1] >= -20) & (v[1:] < -20))[0]
    return (t[ind] * (-v[ind + 1] - 20) + t[ind + 1] * (20 + v[ind])) / (v[ind] - v[ind + 1])


panels = []
for g_syn, t_final in [(0.170, 800.), (0.15, 1600.)]:
    v = simulate(g_syn, t_final)
    t_spikes = spike_times(v, t_final)
    delta = t_spikes - np.floor(t_spikes / T) * T
    panels.append((g_syn, t_final, v, delta))

if __name__ == "__main__":

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))

    for col, (g_syn, t_final, v, delta) in enumerate(panels):
        m_steps = round(t_final / dt)
        t = np.arange(m_steps + 1) * dt
        period_lines = np.arange(1, round(t_final / T) + 2) * T

        ax[0, col].plot(t, v, '-k', linewidth=2)
        for tt in period_lines:
            ax[0, col].axvline(tt, color='r', linestyle=':', linewidth=1)
        ax[0, col].set_xlabel('$t$ [ms]')
        ax[0, col].set_xlim(0, t_final)
        ax[0, col].set_ylim(-100, 50)
        ax[0, col].set_title(rf'$\overline{{g}}_{{\rm syn}}={g_syn}$')

        ax[1, col].plot(np.arange(1, len(delta) + 1), delta, '.k', markersize=10)
        ax[1, col].set_xlabel('spike #')
        ax[1, col].set_xlim(0, len(delta) + 1)
        ax[1, col].set_ylim(0, delta.max() + 1)

    ax[0, 0].set_ylabel('$v$ [mV]')
    ax[1, 0].set_ylabel(r'$\delta$ [ms]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
