import numpy as np
from numpy import exp, tanh
import matplotlib.pyplot as plt

# MATLAB's rng('default'); rng(63806) cannot be bit-reproduced by NumPy,
# so we use our own seed and verify statistically/visually instead of
# expecting an exact match to MATLAB's spike times.
rng = np.random.default_rng(63806)

num_i = 100
v_rev_i = -75.
tau_r_i, tau_peak_i, tau_d_i = 0.5, 0.5, 9.
t_final = 500.
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


tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)


def wb_init_population(i_ext, phi_vec):
    '''vectorized wb_init over a population: each of len(i_ext) WB
    neurons is integrated (Heun/midpoint) independently until its 3rd
    spike, then (v,h,n) is interpolated at phase phi_vec[i] between the
    2nd and 3rd spikes. Faithfully reproduces the matlab source's bug:
    m_tmp is computed from the pre-half-step v, not from v_tmp.'''
    num = len(i_ext)
    max_spikes = 3
    t_final_init = 2000.
    dt_ = 0.01
    dt05_ = dt_ / 2

    v = -70. * np.ones(num)
    m = m_i_inf(v)
    h = h_i_inf(v)
    n = n_i_inf(v)
    t = 0.

    num_spikes = np.zeros(num, dtype=int)
    done = np.zeros(num, dtype=bool)
    t_spikes = np.zeros((num, max_spikes))
    out = np.zeros((num, 3))

    c, g_k, g_na, g_l = 1., 9., 35., 0.1
    v_k, v_na, v_l = -90., 55., -65.

    while np.sum(done) < num and t < t_final_init:
        v_old, h_old, n_old, t_old = v, h, n, t

        v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v) + i_ext) / c
        h_inc = (h_i_inf(v) - h) / tau_h_i(v)
        n_inc = (n_i_inf(v) - n) / tau_n_i(v)

        v_tmp = v + dt05_ * v_inc
        m_tmp = m_i_inf(v)  # faithful port of matlab's bug (uses v, not v_tmp)
        h_tmp = h + dt05_ * h_inc
        n_tmp = n + dt05_ * n_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = (h_i_inf(v_tmp) - h_tmp) / tau_h_i(v_tmp)
        n_inc = (n_i_inf(v_tmp) - n_tmp) / tau_n_i(v_tmp)

        v = v + dt_ * v_inc
        m = m_i_inf(v)
        h = h + dt_ * h_inc
        n = n + dt_ * n_inc
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
        done[ind] = True

    ind = np.where(~done)[0]
    out[ind, 0] = v[ind]
    out[ind, 1] = h[ind]
    out[ind, 2] = n[ind]
    return out


def make_g_ii_fixed_indegree(g_hat_ii, p_ii):
    '''each column (post-synaptic I-cell) keeps exactly round(p_ii*num_i)
    presynaptic connections at the full weight g_hat_ii/(num_i*p_ii),
    rather than each entry independently being present with probability
    p_ii (as in the Erdos-Renyi make_g_ii_random above).'''
    g_ii = g_hat_ii * np.ones((num_i, num_i)) / (num_i * p_ii)
    omit = round(num_i - p_ii * num_i)
    for j in range(num_i):
        drop = rng.choice(num_i, size=omit, replace=False)
        g_ii[drop, j] = 0.
    return g_ii


def simulate(sigma_i, g_hat_ii, p_ii, g_hat_gap, p_gap, t_final=t_final, fixed_indegree=False):
    i_ext_i = 1.5 * (1 + rng.standard_normal(num_i) * sigma_i)

    if fixed_indegree:
        g_ii = make_g_ii_fixed_indegree(g_hat_ii, p_ii)
    else:
        u_ii = rng.random((num_i, num_i))
        g_ii = g_hat_ii * (u_ii < p_ii) / (num_i * p_ii)

    G_gap = np.zeros((num_i, num_i))
    for i in range(num_i - 1):
        for j in range(i + 1, num_i):
            u = rng.random()
            val = (u < p_gap) * g_hat_gap / (p_gap * (num_i - 1))
            G_gap[i, j] = val
            G_gap[j, i] = val
    c = G_gap.sum(axis=1)

    m_steps_ = round(t_final / dt)

    iv = wb_init_population(i_ext_i, rng.random(num_i))
    v_i, h_i, n_i = iv[:, 0], iv[:, 1], iv[:, 2]
    m_i = m_i_inf(v_i)
    q_i, s_i = np.zeros(num_i), np.zeros(num_i)

    t_i_spikes, i_i_spikes = [], []

    for k in range(1, m_steps_ + 1):
        v_i_inc = (0.1 * (-65 - v_i) + 9 * n_i ** 4 * (-90 - v_i) + 35 * m_i ** 3 * h_i * (55 - v_i)
                    + (g_ii.T @ s_i) * (v_rev_i - v_i) + i_ext_i + G_gap @ v_i - c * v_i)
        n_i_inc = (n_i_inf(v_i) - n_i) / tau_n_i(v_i)
        h_i_inc = (h_i_inf(v_i) - h_i) / tau_h_i(v_i)
        q_i_inc = (1 + tanh(v_i / 10)) / 2 * (1 - q_i) / 0.1 - q_i / tau_dq_i
        s_i_inc = q_i * (1 - s_i) / tau_r_i - s_i / tau_d_i

        v_i_tmp = v_i + dt05 * v_i_inc
        n_i_tmp = n_i + dt05 * n_i_inc
        m_i_tmp = m_i_inf(v_i_tmp)
        h_i_tmp = h_i + dt05 * h_i_inc
        q_i_tmp = q_i + dt05 * q_i_inc
        s_i_tmp = s_i + dt05 * s_i_inc

        v_i_inc = (0.1 * (-65 - v_i_tmp) + 9 * n_i_tmp ** 4 * (-90 - v_i_tmp)
                    + 35 * m_i_tmp ** 3 * h_i_tmp * (55 - v_i_tmp)
                    + (g_ii.T @ s_i_tmp) * (v_rev_i - v_i_tmp) + i_ext_i + G_gap @ v_i_tmp - c * v_i_tmp)
        n_i_inc = (n_i_inf(v_i_tmp) - n_i_tmp) / tau_n_i(v_i_tmp)
        h_i_inc = (h_i_inf(v_i_tmp) - h_i_tmp) / tau_h_i(v_i_tmp)
        q_i_inc = (1 + tanh(v_i_tmp / 10)) / 2 * (1 - q_i_tmp) / 0.1 - q_i_tmp / tau_dq_i
        s_i_inc = q_i_tmp * (1 - s_i_tmp) / tau_r_i - s_i_tmp / tau_d_i

        v_i_old = v_i

        v_i = v_i + dt * v_i_inc
        m_i = m_i_inf(v_i)
        h_i = h_i + dt * h_i_inc
        n_i = n_i + dt * n_i_inc
        q_i = q_i + dt * q_i_inc
        s_i = s_i + dt * s_i_inc

        which_i = np.where((v_i_old > -20) & (v_i <= -20))[0]
        if len(which_i) > 0:
            i_i_spikes.extend(which_i.tolist())
            t_i_spikes.extend((((-20 - v_i[which_i]) * (k - 1) * dt + (v_i_old[which_i] + 20) * k * dt)
                                / (-v_i[which_i] + v_i_old[which_i])).tolist())

    return np.array(t_i_spikes), np.array(i_i_spikes)


t_i_spikes, i_i_spikes = simulate(sigma_i=0.03, g_hat_ii=0.5, p_ii=1., g_hat_gap=0.0, p_gap=1.)

if __name__ == "__main__":

    plt.figure(figsize=(8, 4))
    if len(t_i_spikes) > 0:
        plt.plot(t_i_spikes, i_i_spikes, '.k', markersize=2)
    plt.yticks([1, num_i])
    plt.xlabel('$t$ [ms]')
    plt.axis([0, t_final, 0, num_i + 1])
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
