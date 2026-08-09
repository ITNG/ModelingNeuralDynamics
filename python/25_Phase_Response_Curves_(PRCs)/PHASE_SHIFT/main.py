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

i_ext = 0.3

t_final = 300.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

tau_r = 0.5
tau_peak = 0.5
tau_d = 2.
g_syn = 0.1


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


def rtm_init(i_ext, phi0):
    '''find the RTM limit cycle at i_ext (plain-float Heun integration
    until the 5th spike), then interpolate (v, h, n) at phase phi0
    (fraction of the last full period, measured from the 4th spike).'''
    t_final_init = 5000.
    dt_init = 0.001
    dt05_init = dt_init / 2

    v = [-70.]
    m = [m_inf(v[0])]
    h = [h_inf(v[0])]
    n = [n_inf(v[0])]
    t_spikes = []

    k = 0
    t = 0.
    while len(t_spikes) < 5 and t < t_final_init:
        vk, mk, hk, nk = v[k], m[k], h[k], n[k]
        v_inc = (g_k * nk ** 4 * (v_k - vk) + g_na * mk ** 3 * hk * (v_na - vk) + g_l * (v_l - vk) + i_ext) / c
        h_inc = alpha_h(vk) * (1 - hk) - beta_h(vk) * hk
        n_inc = alpha_n(vk) * (1 - nk) - beta_n(vk) * nk

        v_tmp = vk + dt05_init * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = hk + dt05_init * h_inc
        n_tmp = nk + dt05_init * n_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v.append(vk + dt_init * v_inc)
        h.append(hk + dt_init * h_inc)
        n.append(nk + dt_init * n_inc)
        m.append(m_inf(v[-1]))

        if vk >= -20 and v[-1] < -20:
            t_spike = ((k) * dt_init * (-20 - v[-1]) + (k + 1) * dt_init * (20 + vk)) / (vk - v[-1])
            t_spikes.append(t_spike)
        t = (k + 1) * dt_init
        k += 1

    T = t_spikes[4] - t_spikes[3]
    t0 = phi0 * T + t_spikes[3]
    kk = int(t0 / dt_init)
    frac_hi = (t0 - kk * dt_init) / dt_init
    frac_lo = ((kk + 1) * dt_init - t0) / dt_init
    v0 = v[kk + 1] * frac_hi + v[kk] * frac_lo
    h0 = h[kk + 1] * frac_hi + h[kk] * frac_lo
    n0 = n[kk + 1] * frac_hi + n[kk] * frac_lo
    return v0, h0, n0


def simulate_baseline(v0, h0, n0):
    '''free-running RTM neuron, no synaptic input at all.'''
    v = np.zeros(m_steps + 1)
    m = np.zeros(m_steps + 1)
    h = np.zeros(m_steps + 1)
    n = np.zeros(m_steps + 1)
    v[0], h[0], n[0] = v0, h0, n0
    m[0] = m_inf(v0)

    for k in range(m_steps):
        v_inc = (g_k * n[k] ** 4 * (v_k - v[k]) + g_na * m[k] ** 3 * h[k] * (v_na - v[k])
                  + g_l * (v_l - v[k]) + i_ext) / c
        n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
        h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]

        v_tmp = v[k] + dt05 * v_inc
        h_tmp = h[k] + dt05 * h_inc
        n_tmp = n[k] + dt05 * n_inc
        m_tmp = m_inf(v_tmp)

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v[k + 1] = v[k] + dt * v_inc
        m[k + 1] = m_inf(v[k + 1])
        h[k + 1] = h[k] + dt * h_inc
        n[k + 1] = n[k] + dt * n_inc

    return v


def simulate_perturbed(v0, h0, n0, t_pulse):
    '''RTM neuron with a single synaptic input pulse (q set to 1) delivered
    at time t_pulse.'''
    v = np.zeros(m_steps + 1)
    m = np.zeros(m_steps + 1)
    h = np.zeros(m_steps + 1)
    n = np.zeros(m_steps + 1)
    q = np.zeros(m_steps + 1)
    s = np.zeros(m_steps + 1)
    v[0], h[0], n[0] = v0, h0, n0
    m[0] = m_inf(v0)

    for k in range(m_steps):
        v_inc = (g_k * n[k] ** 4 * (v_k - v[k]) + g_na * m[k] ** 3 * h[k] * (v_na - v[k])
                  + g_l * (v_l - v[k]) + g_syn * s[k] * (-v[k]) + i_ext) / c
        n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
        h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]
        q_inc = -q[k] / tau_dq
        s_inc = q[k] * (1 - s[k]) / tau_r - s[k] / tau_d

        v_tmp = v[k] + dt05 * v_inc
        h_tmp = h[k] + dt05 * h_inc
        n_tmp = n[k] + dt05 * n_inc
        m_tmp = m_inf(v_tmp)
        q_tmp = q[k] + dt05 * q_inc
        s_tmp = s[k] + dt05 * s_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + g_syn * s_tmp * (-v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
        q_inc = -q_tmp / tau_dq
        s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d

        v[k + 1] = v[k] + dt * v_inc
        h[k + 1] = h[k] + dt * h_inc
        n[k + 1] = n[k] + dt * n_inc
        m[k + 1] = m_inf(v[k + 1])
        q[k + 1] = q[k] + dt * q_inc
        s[k + 1] = s[k] + dt * s_inc

        if abs((k + 1) * dt - t_pulse) < 1e-6:
            q[k + 1] = 1.

    return v


t = np.arange(m_steps + 1) * dt

v0_1, h0_1, n0_1 = rtm_init(i_ext, 0.2)
v_baseline_1 = simulate_baseline(v0_1, h0_1, n0_1)
v0_2, h0_2, n0_2 = rtm_init(i_ext, 0.2)
v_perturbed_1 = simulate_perturbed(v0_2, h0_2, n0_2, t_pulse=100.)

v0_3, h0_3, n0_3 = rtm_init(i_ext, 0.2)
v_baseline_2 = simulate_baseline(v0_3, h0_3, n0_3)
v0_4, h0_4, n0_4 = rtm_init(i_ext, 0.2)
v_perturbed_2 = simulate_perturbed(v0_4, h0_4, n0_4, t_pulse=120.)

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    axes[0].plot(t, v_baseline_1, '-b', linewidth=4)
    axes[0].axvline(100, linestyle='--', color='k', linewidth=1)
    axes[0].plot(t, v_perturbed_1, '-r', linewidth=1)
    axes[0].set_ylabel('$v$ [mV]')

    axes[1].plot(t, v_baseline_2, '-b', linewidth=4)
    axes[1].axvline(120, linestyle='--', color='k', linewidth=1)
    axes[1].plot(t, v_perturbed_2, '-r', linewidth=1)
    axes[1].set_xlabel('$t$ [ms]')
    axes[1].set_ylabel('$v$ [mV]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
