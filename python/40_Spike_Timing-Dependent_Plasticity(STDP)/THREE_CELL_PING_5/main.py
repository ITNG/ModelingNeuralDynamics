import numpy as np
from numpy import exp, tanh
import matplotlib.pyplot as plt

num_e = 2
num_i = 1
i_ext_e = np.array([0.4, 0.8])
i_ext_i = np.array([0.0])

g_ei = np.array([[0.125], [0.125]])
g_ie = np.array([[0.25, 0.25]])
g_ii = np.array([[0.25]])

v_rev_e, v_rev_i = 0., -75.
tau_r_e, tau_peak_e, tau_d_e = 0.5, 0.5, 3.
tau_r_i, tau_peak_i, tau_d_i = 0.5, 0.5, 9.
t_final = 500.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

# STDP parameters

g_ee = np.zeros((num_e, num_e))
g_ee[0, 1] = 0.05  # initial strengths of the two recurrent E-to-E
g_ee[1, 0] = 0.05  # connections
C = 1.45  # constant in the approximate delta function C*(1+tanh(v_e/10))
K_plus = g_ee.copy()  # max increase in E-to-E strength per spike pair
K_minus = g_ee * 2 / 3  # max decrease in E-to-E strength per spike pair
tau_plus = 10.
tau_minus = 10.
B = 8 * g_ee  # upper bounds on E-to-E synaptic strengths (g_ee is
              # time-dependent, B is fixed)
delta_smooth = np.maximum(g_ee / 2, 1e-6)  # smoothing parameter for the
                                            # soft min/max clamps below


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


def g_ee_derivative(v_e, g_ee_now, a):
    '''STDP increment for the E-to-E weight matrix: a pre-triggered
    potentiation term (y) and a post-triggered depression term (z), each
    softly clamped to [0, B] via a smoothed min/max, and gated by the
    approximate delta functions C*(1+tanh(v_e/10)) of both the
    presynaptic and postsynaptic E-cell.'''
    trace_plus = np.exp(-a / tau_plus) * (1 - np.exp(-5 * a / tau_plus))
    y = g_ee_now + np.diag(trace_plus) @ K_plus
    y = np.minimum(B, y) - delta_smooth / 2 * np.log(1 + np.exp(-2 * np.abs(B - y) / delta_smooth))
    y = y - g_ee_now

    trace_minus = np.exp(-a / tau_minus) * (1 - np.exp(-5 * a / tau_minus))
    z = g_ee_now - K_minus @ np.diag(trace_minus)
    z = np.maximum(0, z) + delta_smooth / 2 * np.log(1 + np.exp(-2 * np.abs(z) / delta_smooth))
    z = z - g_ee_now

    pre = np.tile(tanh(v_e / 10), (num_e, 1))  # pre[i,j] = tanh(v_e[j]/10)
    post = np.tile(tanh(v_e / 10)[:, None], (1, num_e))  # post[i,j] = tanh(v_e[i]/10)
    return C * (1 + pre) * y + C * (1 + post) * z


v_e = np.array([-70., -70.])
m_e, h_e, n_e = m_e_inf(v_e), h_e_inf(v_e), n_e_inf(v_e)
q_e, s_e = np.zeros(num_e), np.zeros(num_e)
a = 100. * np.ones(num_e)

v_i = -75. * np.ones(num_i)
m_i, h_i, n_i = m_i_inf(v_i), h_i_inf(v_i), n_i_inf(v_i)
q_i, s_i = np.zeros(num_i), np.zeros(num_i)

t_e_spikes, i_e_spikes = [], []
t_i_spikes, i_i_spikes = [], []
g_12 = np.zeros(m_steps)
g_21 = np.zeros(m_steps)

for k in range(1, m_steps + 1):
    v_e_inc = (0.1 * (-67 - v_e) + 80 * n_e ** 4 * (-100 - v_e) + 100 * m_e ** 3 * h_e * (50 - v_e)
                + (g_ee.T @ s_e) * (v_rev_e - v_e) + (g_ie.T @ s_i) * (v_rev_i - v_e) + i_ext_e)
    n_e_inc = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)
    h_e_inc = (h_e_inf(v_e) - h_e) / tau_h_e(v_e)
    q_e_inc = (1 + tanh(v_e / 10)) / 2 * (1 - q_e) / 0.1 - q_e / tau_dq_e
    s_e_inc = q_e * (1 - s_e) / tau_r_e - s_e / tau_d_e
    a_inc = 1 - 5 * a * (1 + tanh(v_e / 10))
    g_ee_inc = g_ee_derivative(v_e, g_ee, a)

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
    a_tmp = a + dt05 * a_inc
    g_ee_tmp = g_ee + dt05 * g_ee_inc

    v_i_tmp = v_i + dt05 * v_i_inc
    n_i_tmp = n_i + dt05 * n_i_inc
    m_i_tmp = m_i_inf(v_i_tmp)
    h_i_tmp = h_i + dt05 * h_i_inc
    q_i_tmp = q_i + dt05 * q_i_inc
    s_i_tmp = s_i + dt05 * s_i_inc

    v_e_inc = (0.1 * (-67 - v_e_tmp) + 80 * n_e_tmp ** 4 * (-100 - v_e_tmp)
                + 100 * m_e_tmp ** 3 * h_e_tmp * (50 - v_e_tmp)
                + (g_ee.T @ s_e_tmp) * (v_rev_e - v_e_tmp) + (g_ie.T @ s_i_tmp) * (v_rev_i - v_e_tmp) + i_ext_e)
    n_e_inc = (n_e_inf(v_e_tmp) - n_e_tmp) / tau_n_e(v_e_tmp)
    h_e_inc = (h_e_inf(v_e_tmp) - h_e_tmp) / tau_h_e(v_e_tmp)
    q_e_inc = (1 + tanh(v_e_tmp / 10)) / 2 * (1 - q_e_tmp) / 0.1 - q_e_tmp / tau_dq_e
    s_e_inc = q_e_tmp * (1 - s_e_tmp) / tau_r_e - s_e_tmp / tau_d_e
    a_inc = 1 - 5 * a_tmp * (1 + tanh(v_e_tmp / 10))
    g_ee_inc = g_ee_derivative(v_e_tmp, g_ee_tmp, a_tmp)

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
    a = a + dt * a_inc
    g_ee = g_ee + dt * g_ee_inc

    v_i = v_i + dt * v_i_inc
    m_i = m_i_inf(v_i)
    h_i = h_i + dt * h_i_inc
    n_i = n_i + dt * n_i_inc
    q_i = q_i + dt * q_i_inc
    s_i = s_i + dt * s_i_inc

    which_e = np.where((v_e_old < -40) & (v_e >= -40))[0]
    which_i = np.where((v_i_old < -40) & (v_i >= -40))[0]
    if len(which_e) > 0:
        i_e_spikes.extend(which_e.tolist())
        t_e_spikes.extend((((v_e[which_e] + 40) * (k - 1) * dt + (-v_e_old[which_e] - 40) * k * dt)
                            / (v_e[which_e] - v_e_old[which_e])).tolist())
    if len(which_i) > 0:
        i_i_spikes.extend(which_i.tolist())
        t_i_spikes.extend((((v_i[which_i] + 40) * (k - 1) * dt + (-v_i_old[which_i] - 40) * k * dt)
                            / (v_i[which_i] - v_i_old[which_i])).tolist())

    g_12[k - 1] = g_ee[0, 1]
    g_21[k - 1] = g_ee[1, 0]

t_e_spikes, i_e_spikes = np.array(t_e_spikes), np.array(i_e_spikes)
t_i_spikes, i_i_spikes = np.array(t_i_spikes), np.array(i_i_spikes)

# lag between each spike of E-cell 2 (drive 0.8, index 1) and the next
# following spike of E-cell 1 (drive 0.4, index 0)
t1 = t_e_spikes[i_e_spikes == 0]
t2 = t_e_spikes[i_e_spikes == 1]
lags = []
lag_times = []
for tj in t2:
    later = t1[t1 > tj]
    if len(later) > 0:
        lags.append(later.min() - tj)
        lag_times.append(tj)
lags = np.array(lags)
lag_times = np.array(lag_times)

if __name__ == "__main__":
    t = np.arange(1, m_steps + 1) * dt

    fig, axes = plt.subplots(4, 1, figsize=(8, 9))

    ax = axes[0]
    if len(t_i_spikes) > 0:
        ax.plot(t_i_spikes, i_i_spikes + 1, '.b', markersize=12)
    if len(t_e_spikes) > 0:
        ax.plot(t_e_spikes, i_e_spikes + num_i + 1, '.r', markersize=12)
    ax.plot([0, t_final], [num_i + 0.5, num_i + 0.5], '--k', linewidth=1)
    ax.set_yticks([])
    ax.axis([0, t_final, 0, num_e + num_i + 1])

    axes[1].plot(t, g_12, '-k', linewidth=2)
    axes[1].axis([0, t_final, 0, B[0, 1]])
    axes[1].set_ylabel(r'$\overline{g}_{syn,EE,12}$')

    axes[2].plot(t, g_21, '-k', linewidth=2)
    axes[2].axis([0, t_final, 0, B[1, 0]])
    axes[2].set_ylabel(r'$\overline{g}_{syn,EE,21}$')

    axes[3].plot(lag_times, lags, '.k', markersize=12)
    axes[3].axis([0, t_final, 0, 12])
    axes[3].set_xlabel('$t$ [ms]')
    axes[3].set_ylabel('lag')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
