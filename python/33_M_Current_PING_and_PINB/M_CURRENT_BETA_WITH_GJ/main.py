import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

# MATLAB's rng('default'); rng(63806) cannot be bit-reproduced by NumPy,
# so we use our own seed and verify statistically/visually instead of
# expecting an exact match to MATLAB's spike times.
rng = np.random.default_rng(63806)

num_e = 200
sigma_e = 0.05
i_ext_e = 3.0 * np.ones(num_e) * (1 + sigma_e * rng.standard_normal(num_e))
g_m = 1.0  # M-current conductance; the only parameter of the neuronal
           # model, other than external drive, not explicitly written
           # into the code.
g_hat_gap = 0.2
p_gap = 0.1

t_final = 500.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


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


g_gap = np.zeros((num_e, num_e))
for i in range(num_e - 1):
    u = rng.random(num_e - 1 - i)
    connected = np.where(u < p_gap)[0] + i + 1
    g_gap[i, connected] = g_hat_gap / (num_e - 1) / p_gap
    g_gap[connected, i] = g_gap[i, connected]
g_gap_column_sum = g_gap.sum(axis=0)

iv = rtm_init_with_m_current_population(i_ext_e, rng.random(num_e), g_m)
v_e, h_e, n_e, w = iv[:, 0], iv[:, 1], iv[:, 2], iv[:, 3]
m_e = m_e_inf(v_e)

t_e_spikes, i_e_spikes = [], []

for k in range(1, m_steps + 1):
    v_e_inc = (0.1 * (-67 - v_e) + 80 * n_e ** 4 * (-100 - v_e) + 100 * m_e ** 3 * h_e * (50 - v_e)
                + g_gap @ v_e - g_gap_column_sum * v_e + g_m * w * (-100 - v_e) + i_ext_e)
    n_e_inc = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)
    h_e_inc = (h_e_inf(v_e) - h_e) / tau_h_e(v_e)
    w_inc = (w_inf(v_e) - w) / tau_w(v_e)

    v_e_tmp = v_e + dt05 * v_e_inc
    n_e_tmp = n_e + dt05 * n_e_inc
    m_e_tmp = m_e_inf(v_e_tmp)
    h_e_tmp = h_e + dt05 * h_e_inc
    w_tmp = w + dt05 * w_inc

    v_e_inc = (0.1 * (-67 - v_e_tmp) + 80 * n_e_tmp ** 4 * (-100 - v_e_tmp)
                + 100 * m_e_tmp ** 3 * h_e_tmp * (50 - v_e_tmp)
                + g_gap @ v_e_tmp - g_gap_column_sum * v_e_tmp
                + g_m * w_tmp * (-100 - v_e_tmp) + i_ext_e)
    n_e_inc = (n_e_inf(v_e_tmp) - n_e_tmp) / tau_n_e(v_e_tmp)
    h_e_inc = (h_e_inf(v_e_tmp) - h_e_tmp) / tau_h_e(v_e_tmp)
    # note: matlab's source does not recompute w_inc from the tmp state
    # here, so the corrector step below faithfully reuses the predictor
    # w_inc (an Euler, not midpoint, update for w -- a quirk of the
    # original code).

    v_e_old = v_e

    v_e = v_e + dt * v_e_inc
    m_e = m_e_inf(v_e)
    h_e = h_e + dt * h_e_inc
    n_e = n_e + dt * n_e_inc
    w = w + dt * w_inc

    which_e = np.where((v_e_old > -20) & (v_e <= -20))[0]
    if len(which_e) > 0:
        i_e_spikes.extend(which_e.tolist())
        t_e_spikes.extend((((-20 - v_e[which_e]) * (k - 1) * dt + (v_e_old[which_e] + 20) * k * dt)
                            / (-v_e[which_e] + v_e_old[which_e])).tolist())

t_e_spikes, i_e_spikes = np.array(t_e_spikes), np.array(i_e_spikes)

if __name__ == "__main__":
    f_hat_e = round(len(t_e_spikes) / num_e / t_final * 1000)
    print(f"f_hat_e = {f_hat_e}")

    fig, ax = plt.subplots(figsize=(8, 4))
    if len(t_e_spikes) > 0:
        ax.plot(t_e_spikes, i_e_spikes, '.r', markersize=2)
    ax.axis([0, t_final, 0, num_e + 1])
    ax.set_xlabel('$t$ [ms]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
