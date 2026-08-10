import numpy as np
from numpy import exp, tanh
import matplotlib.pyplot as plt

v_rev_i = -75.
tau_r_i, tau_peak_i = 0.5, 0.5

t_final = 200.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


def m_i_inf(v):
    alpha_m = 0.1 * (v + 35) / (1 - exp(-(v + 35) / 10))
    beta_m = 4. * exp(-(v + 60) / 18)
    return alpha_m / (alpha_m + beta_m)


def h_i_inf(v):
    alpha_h = 0.07 * exp(-(v + 58) / 20)
    beta_h = 1. / (exp(-0.1 * (v + 28)) + 1)
    return alpha_h / (alpha_h + beta_h)


def tau_h_i(v):
    alpha_h = 0.07 * exp(-(v + 58) / 20)
    beta_h = 1. / (exp(-0.1 * (v + 28)) + 1)
    return 1. / (alpha_h + beta_h) / 5


def n_i_inf(v):
    alpha_n = -0.01 * (v + 34) / (exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * exp(-(v + 44) / 80)
    return alpha_n / (alpha_n + beta_n)


def tau_n_i(v):
    alpha_n = -0.01 * (v + 34) / (exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * exp(-(v + 44) / 80)
    return 1. / (alpha_n + beta_n) / 5


def tau_peak_function(tau_d, tau_r, tau_d_q):
    dt_ = 0.01
    dt05_ = dt_ / 2
    s, t = 0., 0.
    s_inc = exp(-t / tau_d_q) * (1 - s) / tau_r - s * tau_d
    while s_inc > 0:
        t_old, s_inc_old = t, s_inc
        s_tmp = s + dt05_ * s_inc
        s_inc_tmp = exp(-(t + dt05_) / tau_d_q) * (1 - s_tmp) / tau_r - s_tmp / tau_d
        s = s + dt_ * s_inc_tmp
        t = t + dt_
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


def simulate(i_ext_i, g_ii, tau_d_i, record_trace=False):
    tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)

    v_i = [-75.]
    h_i, n_i, q_i, s_i = 0.1, 0.1, 0., 0.
    t_i_spikes = []

    for k in range(m_steps):
        t_old, t_new = k * dt, (k + 1) * dt
        vk = v_i[k]

        v_inc = (0.1 * (-65 - vk) + 9 * n_i ** 4 * (-90 - vk) + 35 * m_i_inf(vk) ** 3 * h_i * (55 - vk)
                  + g_ii * s_i * (v_rev_i - vk) + i_ext_i)
        n_inc = (n_i_inf(vk) - n_i) / tau_n_i(vk)
        h_inc = (h_i_inf(vk) - h_i) / tau_h_i(vk)
        q_inc = (1 + tanh(vk / 10)) / 2 * (1 - q_i) / 0.1 - q_i / tau_dq_i
        s_inc = q_i * (1 - s_i) / tau_r_i - s_i / tau_d_i

        v_tmp = vk + dt05 * v_inc
        n_tmp = n_i + dt05 * n_inc
        h_tmp = h_i + dt05 * h_inc
        q_tmp = q_i + dt05 * q_inc
        s_tmp = s_i + dt05 * s_inc

        v_inc = (0.1 * (-65 - v_tmp) + 9 * n_tmp ** 4 * (-90 - v_tmp) + 35 * m_i_inf(v_tmp) ** 3 * h_tmp * (55 - v_tmp)
                  + g_ii * s_tmp * (v_rev_i - v_tmp) + i_ext_i)
        n_inc = (n_i_inf(v_tmp) - n_tmp) / tau_n_i(v_tmp)
        h_inc = (h_i_inf(v_tmp) - h_tmp) / tau_h_i(v_tmp)
        q_inc = (1 + tanh(v_tmp / 10)) / 2 * (1 - q_tmp) / 0.1 - q_tmp / tau_dq_i
        s_inc = q_tmp * (1 - s_tmp) / tau_r_i - s_tmp / tau_d_i

        v_i.append(vk + dt * v_inc)
        h_i = h_i + dt * h_inc
        n_i = n_i + dt * n_inc
        q_i = q_i + dt * q_inc
        s_i = s_i + dt * s_inc

        if v_i[k + 1] < -20 and vk >= -20:
            tt = (t_old * (-20 - v_i[k + 1]) + t_new * (vk + 20)) / (vk - v_i[k + 1])
            t_i_spikes.append(tt)

    period = t_i_spikes[-1] - t_i_spikes[-2]
    v_i = np.array(v_i) if record_trace else None
    return period, v_i


base_period, v_trace = simulate(i_ext_i=1.5, g_ii=0.5, tau_d_i=9., record_trace=True)

period_reduced_I, _ = simulate(i_ext_i=1.5 * 0.99, g_ii=0.5, tau_d_i=9.)
pct_i_ext = (base_period - period_reduced_I) / base_period * 100

period_raised_g_ii, _ = simulate(i_ext_i=1.5, g_ii=0.5 * 1.01, tau_d_i=9.)
pct_g_ii = (base_period - period_raised_g_ii) / base_period * 100

period_raised_tau_d, _ = simulate(i_ext_i=1.5, g_ii=0.5, tau_d_i=9. * 1.01)
pct_tau_d = (base_period - period_raised_tau_d) / base_period * 100

if __name__ == "__main__":
    print("base_period", base_period)
    print("pct_i_ext", pct_i_ext)
    print("pct_g_ii", pct_g_ii)
    print("pct_tau_d", pct_tau_d)

    plt.figure(figsize=(8, 4))
    t = np.arange(m_steps + 1) * dt
    plt.plot(t, v_trace, '-b', linewidth=2)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('$v$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
