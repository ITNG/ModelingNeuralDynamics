import numpy as np
from numpy import exp, tanh
import matplotlib.pyplot as plt

# MATLAB's rng('default'); rng(63806) cannot be bit-reproduced by NumPy,
# so we use our own seed and verify statistically/visually instead of
# expecting an exact match to MATLAB's spike times.
rng = np.random.default_rng(63806)

num_e = 200
num_i = 50
sigma_e = 0.05
i_ext_e = 3.0 * np.ones(num_e) * (1 + sigma_e * rng.standard_normal(num_e))
sigma_i = 0.05
i_ext_i = 0.7 * np.ones(num_i) * (1 + sigma_i * rng.standard_normal(num_i))
g_m = 1.0  # M-current conductance; the only parameter of the neuronal
           # model, other than external drive, not explicitly written
           # into the code.
g_hat_ee, g_hat_ei, g_hat_ie, g_hat_ii = 0., 0.5, 0.5, 0.5
p_ee, p_ei, p_ie, p_ii = 0.5, 0.5, 0.5, 0.5

v_rev_e, v_rev_i = 0., -75.
tau_r_e, tau_peak_e, tau_d_e = 0.5, 0.5, 3.
tau_r_i, tau_peak_i, tau_d_i = 0.5, 0.5, 9.
t_final = 1000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


# ------------------------------------------------------------- E cell (RTM)


def m_e_inf(v):
    alpha_m = 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))
    beta_m = 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)
    return alpha_m / (alpha_m + beta_m)


def h_e_inf(v):
    alpha_h = 0.128 * exp(-(v + 50) / 18)
    beta_h = 4. / (1 + exp(-(v + 27) / 5))
    return alpha_h / (alpha_h + beta_h)


def tau_h_e(v):
    alpha_h = 0.128 * exp(-(v + 50) / 18)
    beta_h = 4. / (1 + exp(-(v + 27) / 5))
    return 1. / (alpha_h + beta_h)


def n_e_inf(v):
    alpha_n = 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))
    beta_n = 0.5 * exp(-(v + 57) / 40)
    return alpha_n / (alpha_n + beta_n)


def tau_n_e(v):
    alpha_n = 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))
    beta_n = 0.5 * exp(-(v + 57) / 40)
    return 1. / (alpha_n + beta_n)


def w_inf(v):
    return 1. / (1 + exp(-(v + 35) / 10))


def tau_w(v):
    return 400. / (3.3 * exp((v + 35) / 20) + exp(-(v + 35) / 20))


# ------------------------------------------------------------- I cell (WB)


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


# ------------------------------------------------------------------- shared


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


tau_dq_e = tau_d_q_function(tau_d_e, tau_r_e, tau_peak_e)
tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)


def rtm_init_with_m_current_population(i_ext, phi_vec, g_m):
    '''vectorized RTM-with-M-current init over a population: each of
    len(i_ext) neurons is integrated (Heun/midpoint) independently until
    its 3rd spike, then (v,h,n,w) is interpolated at phase phi_vec[i]
    between the 2nd and 3rd spikes. Faithfully reproduces the matlab
    source's bug: m_tmp is computed from the pre-half-step v, not v_tmp.'''
    num = len(i_ext)
    max_spikes = 3
    t_final_init = 2000.
    dt_ = 0.01
    dt05_ = dt_ / 2

    v = -70. * np.ones(num)
    m = m_e_inf(v)
    h = h_e_inf(v)
    n = n_e_inf(v)
    w = np.zeros(num)
    t = 0.

    num_spikes = np.zeros(num, dtype=int)
    done = np.zeros(num, dtype=bool)
    t_spikes = np.zeros((num, max_spikes))
    out = np.zeros((num, 4))

    c, g_k, g_na, g_l = 1., 80., 100., 0.1
    v_k, v_na, v_l = -100., 50., -67.

    while np.sum(done) < num and t < t_final_init:
        v_old, h_old, n_old, w_old, t_old = v, h, n, w, t

        v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v)
                  + g_l * (v_l - v) + g_m * w * (v_k - v) + i_ext) / c
        h_inc = (h_e_inf(v) - h) / tau_h_e(v)
        n_inc = (n_e_inf(v) - n) / tau_n_e(v)
        w_inc = (w_inf(v) - w) / tau_w(v)

        v_tmp = v + dt05_ * v_inc
        m_tmp = m_e_inf(v)  # faithful port of matlab's bug (uses v, not v_tmp)
        h_tmp = h + dt05_ * h_inc
        n_tmp = n + dt05_ * n_inc
        w_tmp = w + dt05_ * w_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + g_m * w_tmp * (v_k - v_tmp) + i_ext) / c
        h_inc = (h_e_inf(v_tmp) - h_tmp) / tau_h_e(v_tmp)
        n_inc = (n_e_inf(v_tmp) - n_tmp) / tau_n_e(v_tmp)
        w_inc = (w_inf(v_tmp) - w_tmp) / tau_w(v_tmp)

        v = v + dt_ * v_inc
        m = m_e_inf(v)
        h = h + dt_ * h_inc
        n = n + dt_ * n_inc
        w = w + dt_ * w_inc
        t = t + dt_

        ind = np.where((v_old >= -20) & (v < -20))[0]
        for k in ind:
            num_spikes[k] += 1
            if num_spikes[k] <= max_spikes:
                t_spikes[k, num_spikes[k] - 1] = (t_old * (-20 - v[k]) + t * (v_old[k] + 20)) / (v_old[k] - v[k])

        thr = t_spikes[:, max_spikes - 1] + phi_vec * (t_spikes[:, max_spikes - 1] - t_spikes[:, max_spikes - 2])
        ind = np.where((num_spikes == max_spikes) & (t > thr) & (t_old <= thr) & (~done))[0]
        for k in ind:
            out[k, 0] = (v_old[k] * (t - thr[k]) + v[k] * (thr[k] - t_old)) / dt_
            out[k, 1] = (h_old[k] * (t - thr[k]) + h[k] * (thr[k] - t_old)) / dt_
            out[k, 2] = (n_old[k] * (t - thr[k]) + n[k] * (thr[k] - t_old)) / dt_
            out[k, 3] = (w_old[k] * (t - thr[k]) + w[k] * (thr[k] - t_old)) / dt_
        done[ind] = True

    ind = np.where(~done)[0]
    out[ind, 0] = v[ind]
    out[ind, 1] = h[ind]
    out[ind, 2] = n[ind]
    out[ind, 3] = w[ind]
    return out


def make_random_connectivity(p_ee, p_ei, p_ie, p_ii):
    u_ee = rng.random((num_e, num_e))
    u_ei = rng.random((num_e, num_i))
    u_ie = rng.random((num_i, num_e))
    u_ii = rng.random((num_i, num_i))
    g_ee = g_hat_ee * (u_ee < p_ee) / (num_e * p_ee) if p_ee > 0 else np.zeros((num_e, num_e))
    g_ei = g_hat_ei * (u_ei < p_ei) / (num_e * p_ei)
    g_ie = g_hat_ie * (u_ie < p_ie) / (num_i * p_ie)
    g_ii = g_hat_ii * (u_ii < p_ii) / (num_i * p_ii)
    return g_ee, g_ei, g_ie, g_ii


def simulate(g_ee, g_ei, g_ie, g_ii, t_final=t_final, track_cell=3):
    m_steps_ = round(t_final / dt)

    iv = rtm_init_with_m_current_population(i_ext_e, rng.random(num_e), g_m)
    v_e, h_e, n_e, w = iv[:, 0], iv[:, 1], iv[:, 2], iv[:, 3]
    m_e = m_e_inf(v_e)
    q_e, s_e = np.zeros(num_e), np.zeros(num_e)

    v_i = -75. * np.ones(num_i)
    m_i, h_i, n_i = m_i_inf(v_i), h_i_inf(v_i), n_i_inf(v_i)
    q_i, s_i = np.zeros(num_i), np.zeros(num_i)

    t_e_spikes, i_e_spikes = [], []
    t_i_spikes, i_i_spikes = [], []
    v_plot = np.zeros(m_steps_)

    for k in range(1, m_steps_ + 1):
        v_e_inc = (0.1 * (-67 - v_e) + 80 * n_e ** 4 * (-100 - v_e) + 100 * m_e ** 3 * h_e * (50 - v_e)
                    + (g_ee.T @ s_e) * (v_rev_e - v_e) + (g_ie.T @ s_i) * (v_rev_i - v_e)
                    + g_m * w * (-100 - v_e) + i_ext_e)
        n_e_inc = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)
        h_e_inc = (h_e_inf(v_e) - h_e) / tau_h_e(v_e)
        q_e_inc = (1 + tanh(v_e / 10)) / 2 * (1 - q_e) / 0.1 - q_e / tau_dq_e
        s_e_inc = q_e * (1 - s_e) / tau_r_e - s_e / tau_d_e
        w_inc = (w_inf(v_e) - w) / tau_w(v_e)

        v_i_inc = (0.1 * (-65 - v_i) + 9 * n_i ** 4 * (-90 - v_i) + 35 * m_i ** 3 * h_i * (55 - v_i)
                    + (g_ei.T @ s_e) * (v_rev_e - v_i) + (g_ii.T @ s_i) * (v_rev_i - v_i) + i_ext_i)
        n_i_inc = (n_i_inf(v_i) - n_i) / tau_n_i(v_i)
        h_i_inc = (h_i_inf(v_i) - h_i) / tau_h_i(v_i)
        q_i_inc = (1 + tanh(v_i / 10)) / 2 * (1 - q_i) / 0.1 - q_i / tau_dq_i
        s_i_inc = q_i * (1 - s_i) / tau_r_i - s_i / tau_d_i

        v_e_tmp = v_e + dt05 * v_e_inc
        n_e_tmp = n_e + dt05 * n_e_inc
        m_e_tmp = m_e_inf(v_e_tmp)
        h_e_tmp = h_e + dt05 * h_e_inc
        q_e_tmp = q_e + dt05 * q_e_inc
        s_e_tmp = s_e + dt05 * s_e_inc
        w_tmp = w + dt05 * w_inc

        v_i_tmp = v_i + dt05 * v_i_inc
        n_i_tmp = n_i + dt05 * n_i_inc
        m_i_tmp = m_i_inf(v_i_tmp)
        h_i_tmp = h_i + dt05 * h_i_inc
        q_i_tmp = q_i + dt05 * q_i_inc
        s_i_tmp = s_i + dt05 * s_i_inc

        v_e_inc = (0.1 * (-67 - v_e_tmp) + 80 * n_e_tmp ** 4 * (-100 - v_e_tmp)
                    + 100 * m_e_tmp ** 3 * h_e_tmp * (50 - v_e_tmp)
                    + (g_ee.T @ s_e_tmp) * (v_rev_e - v_e_tmp) + (g_ie.T @ s_i_tmp) * (v_rev_i - v_e_tmp)
                    + g_m * w_tmp * (-100 - v_e_tmp) + i_ext_e)
        n_e_inc = (n_e_inf(v_e_tmp) - n_e_tmp) / tau_n_e(v_e_tmp)
        h_e_inc = (h_e_inf(v_e_tmp) - h_e_tmp) / tau_h_e(v_e_tmp)
        q_e_inc = (1 + tanh(v_e_tmp / 10)) / 2 * (1 - q_e_tmp) / 0.1 - q_e_tmp / tau_dq_e
        s_e_inc = q_e_tmp * (1 - s_e_tmp) / tau_r_e - s_e_tmp / tau_d_e
        w_inc = (w_inf(v_e_tmp) - w_tmp) / tau_w(v_e_tmp)

        v_i_inc = (0.1 * (-65 - v_i_tmp) + 9 * n_i_tmp ** 4 * (-90 - v_i_tmp)
                    + 35 * m_i_tmp ** 3 * h_i_tmp * (55 - v_i_tmp)
                    + (g_ei.T @ s_e_tmp) * (v_rev_e - v_i_tmp) + (g_ii.T @ s_i_tmp) * (v_rev_i - v_i_tmp) + i_ext_i)
        n_i_inc = (n_i_inf(v_i_tmp) - n_i_tmp) / tau_n_i(v_i_tmp)
        h_i_inc = (h_i_inf(v_i_tmp) - h_i_tmp) / tau_h_i(v_i_tmp)
        q_i_inc = (1 + tanh(v_i_tmp / 10)) / 2 * (1 - q_i_tmp) / 0.1 - q_i_tmp / tau_dq_i
        s_i_inc = q_i_tmp * (1 - s_i_tmp) / tau_r_i - s_i_tmp / tau_d_i

        v_e_old, v_i_old = v_e, v_i

        v_e = v_e + dt * v_e_inc
        m_e = m_e_inf(v_e)
        h_e = h_e + dt * h_e_inc
        n_e = n_e + dt * n_e_inc
        q_e = q_e + dt * q_e_inc
        s_e = s_e + dt * s_e_inc
        w = w + dt * w_inc

        v_i = v_i + dt * v_i_inc
        m_i = m_i_inf(v_i)
        h_i = h_i + dt * h_i_inc
        n_i = n_i + dt * n_i_inc
        q_i = q_i + dt * q_i_inc
        s_i = s_i + dt * s_i_inc

        which_e = np.where((v_e_old > -20) & (v_e <= -20))[0]
        which_i = np.where((v_i_old > -20) & (v_i <= -20))[0]
        if len(which_e) > 0:
            i_e_spikes.extend(which_e.tolist())
            t_e_spikes.extend((((-20 - v_e[which_e]) * (k - 1) * dt + (v_e_old[which_e] + 20) * k * dt)
                                / (-v_e[which_e] + v_e_old[which_e])).tolist())
        if len(which_i) > 0:
            i_i_spikes.extend(which_i.tolist())
            t_i_spikes.extend((((-20 - v_i[which_i]) * (k - 1) * dt + (v_i_old[which_i] + 20) * k * dt)
                                / (-v_i[which_i] + v_i_old[which_i])).tolist())

        v_plot[k - 1] = v_e[track_cell]

    return (np.array(t_e_spikes), np.array(i_e_spikes),
            np.array(t_i_spikes), np.array(i_i_spikes), v_plot)


def run_smoke(t_final_run=t_final):
    connectivity = make_random_connectivity(p_ee, p_ei, p_ie, p_ii)
    return simulate(*connectivity, t_final=t_final_run)


if __name__ == "__main__":
    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, v_plot = run_smoke()
    f_hat_e = round(len(t_e_spikes) / num_e / t_final * 1000)
    f_hat_i = round(len(t_i_spikes) / num_i / t_final * 1000)
    print(f"f_hat_e = {f_hat_e}")
    print(f"f_hat_i = {f_hat_i}")

    fig, axes = plt.subplots(2, 1, figsize=(8, 6))

    ax = axes[0]
    if len(t_i_spikes) > 0:
        ax.plot(t_i_spikes, i_i_spikes, '.b', markersize=2)
    if len(t_e_spikes) > 0:
        ax.plot(t_e_spikes, i_e_spikes + num_i, '.r', markersize=2)
    ax.plot([0, t_final], [num_i + 0.5, num_i + 0.5], '--k', linewidth=1)
    ax.set_yticks([num_i, num_e + num_i])
    ax.axis([0, t_final, 0, num_e + num_i + 1])

    axes[1].plot(np.arange(1, m_steps + 1) * dt, v_plot, '-k', linewidth=2)
    axes[1].set_xlabel('$t$ [ms]')
    axes[1].set_ylabel('$v$ [mV]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
