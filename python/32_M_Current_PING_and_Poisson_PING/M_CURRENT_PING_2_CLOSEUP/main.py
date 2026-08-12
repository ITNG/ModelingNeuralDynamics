import math

import numpy as np
from numpy import exp, tanh
import matplotlib.pyplot as plt
from numba import njit
from numba.typed import List

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


# ------------------------------------------------------- scalar njit mirrors
# Used only by the compiled main-loop below (numba needs scalar math.exp
# calls, not vectorized numpy, to compile a tight per-neuron loop -- see
# the array-based originals above for the functions these mirror exactly).


@njit
def _m_e_inf_s(v):
    alpha_m = 0.32 * (v + 54) / (1 - math.exp(-(v + 54) / 4))
    beta_m = 0.28 * (v + 27) / (math.exp((v + 27) / 5) - 1)
    return alpha_m / (alpha_m + beta_m)


@njit
def _h_e_inf_s(v):
    alpha_h = 0.128 * math.exp(-(v + 50) / 18)
    beta_h = 4. / (1 + math.exp(-(v + 27) / 5))
    return alpha_h / (alpha_h + beta_h)


@njit
def _tau_h_e_s(v):
    alpha_h = 0.128 * math.exp(-(v + 50) / 18)
    beta_h = 4. / (1 + math.exp(-(v + 27) / 5))
    return 1. / (alpha_h + beta_h)


@njit
def _n_e_inf_s(v):
    alpha_n = 0.032 * (v + 52) / (1 - math.exp(-(v + 52) / 5))
    beta_n = 0.5 * math.exp(-(v + 57) / 40)
    return alpha_n / (alpha_n + beta_n)


@njit
def _tau_n_e_s(v):
    alpha_n = 0.032 * (v + 52) / (1 - math.exp(-(v + 52) / 5))
    beta_n = 0.5 * math.exp(-(v + 57) / 40)
    return 1. / (alpha_n + beta_n)


@njit
def _w_inf_s(v):
    return 1. / (1 + math.exp(-(v + 35) / 10))


@njit
def _tau_w_s(v):
    return 400. / (3.3 * math.exp((v + 35) / 20) + math.exp(-(v + 35) / 20))


@njit
def _m_i_inf_s(v):
    alpha_m = 0.1 * (v + 35) / (1 - math.exp(-(v + 35) / 10))
    beta_m = 4. * math.exp(-(v + 60) / 18)
    return alpha_m / (alpha_m + beta_m)


@njit
def _h_i_inf_s(v):
    alpha_h = 0.07 * math.exp(-(v + 58) / 20)
    beta_h = 1. / (math.exp(-0.1 * (v + 28)) + 1)
    return alpha_h / (alpha_h + beta_h)


@njit
def _tau_h_i_s(v):
    alpha_h = 0.07 * math.exp(-(v + 58) / 20)
    beta_h = 1. / (math.exp(-0.1 * (v + 28)) + 1)
    return 1. / (alpha_h + beta_h) / 5.


@njit
def _n_i_inf_s(v):
    alpha_n = -0.01 * (v + 34) / (math.exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * math.exp(-(v + 44) / 80)
    return alpha_n / (alpha_n + beta_n)


@njit
def _tau_n_i_s(v):
    alpha_n = -0.01 * (v + 34) / (math.exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * math.exp(-(v + 44) / 80)
    return 1. / (alpha_n + beta_n) / 5.


@njit
def _m_current_ping_step_loop(m_steps, dt, dt05, num_e, num_i, g_m,
                               v_rev_e, v_rev_i, tau_r_e, tau_d_e, tau_dq_e,
                               tau_r_i, tau_d_i, tau_dq_i,
                               i_ext_e, i_ext_i, g_ee, g_ei, g_ie, g_ii,
                               v_e, h_e, n_e, m_e, q_e, s_e, w,
                               v_i, h_i, n_i, m_i, q_i, s_i, v_plot, track_cell):
    e_times = List.empty_list(np.float64)
    e_indices = List.empty_list(np.int64)
    i_times = List.empty_list(np.float64)
    i_indices = List.empty_list(np.int64)

    dve = np.empty(num_e); dne = np.empty(num_e); dhe = np.empty(num_e)
    dqe = np.empty(num_e); dse = np.empty(num_e); dw = np.empty(num_e)
    dvi = np.empty(num_i); dni = np.empty(num_i); dhi = np.empty(num_i)
    dqi = np.empty(num_i); dsi = np.empty(num_i)

    ve_m = np.empty(num_e); ne_m = np.empty(num_e); me_m = np.empty(num_e)
    he_m = np.empty(num_e); qe_m = np.empty(num_e); se_m = np.empty(num_e); w_m = np.empty(num_e)
    vi_m = np.empty(num_i); ni_m = np.empty(num_i); mi_m = np.empty(num_i)
    hi_m = np.empty(num_i); qi_m = np.empty(num_i); si_m = np.empty(num_i)

    ve_old = np.empty(num_e)
    vi_old = np.empty(num_i)
    ee_term = np.empty(num_e); ie_term = np.empty(num_e)
    ei_term = np.empty(num_i); ii_term = np.empty(num_i)

    for step in range(m_steps):
        k = step + 1

        for i in range(num_e):
            acc = 0.0
            for j in range(num_e):
                acc += g_ee[j, i] * s_e[j]
            ee_term[i] = acc
            acc = 0.0
            for j in range(num_i):
                acc += g_ie[j, i] * s_i[j]
            ie_term[i] = acc
        for i in range(num_i):
            acc = 0.0
            for j in range(num_e):
                acc += g_ei[j, i] * s_e[j]
            ei_term[i] = acc
            acc = 0.0
            for j in range(num_i):
                acc += g_ii[j, i] * s_i[j]
            ii_term[i] = acc

        for j in range(num_e):
            v = v_e[j]
            dve[j] = (0.1 * (-67 - v) + 80 * n_e[j] ** 4 * (-100 - v)
                      + 100 * m_e[j] ** 3 * h_e[j] * (50 - v)
                      + ee_term[j] * (v_rev_e - v) + ie_term[j] * (v_rev_i - v)
                      + g_m * w[j] * (-100 - v) + i_ext_e[j])
            dne[j] = (_n_e_inf_s(v) - n_e[j]) / _tau_n_e_s(v)
            dhe[j] = (_h_e_inf_s(v) - h_e[j]) / _tau_h_e_s(v)
            th = math.tanh(v / 10)
            dqe[j] = (1 + th) / 2 * (1 - q_e[j]) / 0.1 - q_e[j] / tau_dq_e
            dse[j] = q_e[j] * (1 - s_e[j]) / tau_r_e - s_e[j] / tau_d_e
            dw[j] = (_w_inf_s(v) - w[j]) / _tau_w_s(v)
        for j in range(num_i):
            v = v_i[j]
            dvi[j] = (0.1 * (-65 - v) + 9 * n_i[j] ** 4 * (-90 - v)
                      + 35 * m_i[j] ** 3 * h_i[j] * (55 - v)
                      + ei_term[j] * (v_rev_e - v) + ii_term[j] * (v_rev_i - v) + i_ext_i[j])
            dni[j] = (_n_i_inf_s(v) - n_i[j]) / _tau_n_i_s(v)
            dhi[j] = (_h_i_inf_s(v) - h_i[j]) / _tau_h_i_s(v)
            th = math.tanh(v / 10)
            dqi[j] = (1 + th) / 2 * (1 - q_i[j]) / 0.1 - q_i[j] / tau_dq_i
            dsi[j] = q_i[j] * (1 - s_i[j]) / tau_r_i - s_i[j] / tau_d_i

        for j in range(num_e):
            ve_m[j] = v_e[j] + dt05 * dve[j]
            ne_m[j] = n_e[j] + dt05 * dne[j]
            me_m[j] = _m_e_inf_s(ve_m[j])
            he_m[j] = h_e[j] + dt05 * dhe[j]
            qe_m[j] = q_e[j] + dt05 * dqe[j]
            se_m[j] = s_e[j] + dt05 * dse[j]
            w_m[j] = w[j] + dt05 * dw[j]
        for j in range(num_i):
            vi_m[j] = v_i[j] + dt05 * dvi[j]
            ni_m[j] = n_i[j] + dt05 * dni[j]
            mi_m[j] = _m_i_inf_s(vi_m[j])
            hi_m[j] = h_i[j] + dt05 * dhi[j]
            qi_m[j] = q_i[j] + dt05 * dqi[j]
            si_m[j] = s_i[j] + dt05 * dsi[j]

        for i in range(num_e):
            acc = 0.0
            for j in range(num_e):
                acc += g_ee[j, i] * se_m[j]
            ee_term[i] = acc
            acc = 0.0
            for j in range(num_i):
                acc += g_ie[j, i] * si_m[j]
            ie_term[i] = acc
        for i in range(num_i):
            acc = 0.0
            for j in range(num_e):
                acc += g_ei[j, i] * se_m[j]
            ei_term[i] = acc
            acc = 0.0
            for j in range(num_i):
                acc += g_ii[j, i] * si_m[j]
            ii_term[i] = acc

        for j in range(num_e):
            v = ve_m[j]
            dve[j] = (0.1 * (-67 - v) + 80 * ne_m[j] ** 4 * (-100 - v)
                      + 100 * me_m[j] ** 3 * he_m[j] * (50 - v)
                      + ee_term[j] * (v_rev_e - v) + ie_term[j] * (v_rev_i - v)
                      + g_m * w_m[j] * (-100 - v) + i_ext_e[j])
            dne[j] = (_n_e_inf_s(v) - ne_m[j]) / _tau_n_e_s(v)
            dhe[j] = (_h_e_inf_s(v) - he_m[j]) / _tau_h_e_s(v)
            th = math.tanh(v / 10)
            dqe[j] = (1 + th) / 2 * (1 - qe_m[j]) / 0.1 - qe_m[j] / tau_dq_e
            dse[j] = qe_m[j] * (1 - se_m[j]) / tau_r_e - se_m[j] / tau_d_e
            dw[j] = (_w_inf_s(v) - w_m[j]) / _tau_w_s(v)
        for j in range(num_i):
            v = vi_m[j]
            dvi[j] = (0.1 * (-65 - v) + 9 * ni_m[j] ** 4 * (-90 - v)
                      + 35 * mi_m[j] ** 3 * hi_m[j] * (55 - v)
                      + ei_term[j] * (v_rev_e - v) + ii_term[j] * (v_rev_i - v) + i_ext_i[j])
            dni[j] = (_n_i_inf_s(v) - ni_m[j]) / _tau_n_i_s(v)
            dhi[j] = (_h_i_inf_s(v) - hi_m[j]) / _tau_h_i_s(v)
            th = math.tanh(v / 10)
            dqi[j] = (1 + th) / 2 * (1 - qi_m[j]) / 0.1 - qi_m[j] / tau_dq_i
            dsi[j] = qi_m[j] * (1 - si_m[j]) / tau_r_i - si_m[j] / tau_d_i

        for j in range(num_e):
            ve_old[j] = v_e[j]
        for j in range(num_i):
            vi_old[j] = v_i[j]

        for j in range(num_e):
            v_e[j] = v_e[j] + dt * dve[j]
            m_e[j] = _m_e_inf_s(v_e[j])
            h_e[j] = h_e[j] + dt * dhe[j]
            n_e[j] = n_e[j] + dt * dne[j]
            q_e[j] = q_e[j] + dt * dqe[j]
            s_e[j] = s_e[j] + dt * dse[j]
            w[j] = w[j] + dt * dw[j]
        for j in range(num_i):
            v_i[j] = v_i[j] + dt * dvi[j]
            m_i[j] = _m_i_inf_s(v_i[j])
            h_i[j] = h_i[j] + dt * dhi[j]
            n_i[j] = n_i[j] + dt * dni[j]
            q_i[j] = q_i[j] + dt * dqi[j]
            s_i[j] = s_i[j] + dt * dsi[j]

        for j in range(num_e):
            if ve_old[j] > -20 and v_e[j] <= -20:
                e_indices.append(j)
                e_times.append(((-20 - v_e[j]) * step * dt + (ve_old[j] + 20) * k * dt)
                                / (ve_old[j] - v_e[j]))
        for j in range(num_i):
            if vi_old[j] > -20 and v_i[j] <= -20:
                i_indices.append(j)
                i_times.append(((-20 - v_i[j]) * step * dt + (vi_old[j] + 20) * k * dt)
                                / (vi_old[j] - v_i[j]))

        v_plot[step] = v_e[track_cell]

    return e_times, e_indices, i_times, i_indices


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

    v_plot = np.zeros(m_steps_)

    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes = _m_current_ping_step_loop(
        m_steps_, dt, dt05, num_e, num_i, g_m,
        v_rev_e, v_rev_i, tau_r_e, tau_d_e, tau_dq_e,
        tau_r_i, tau_d_i, tau_dq_i,
        i_ext_e, i_ext_i, g_ee, g_ei, g_ie, g_ii,
        v_e, h_e, n_e, m_e, q_e, s_e, w,
        v_i, h_i, n_i, m_i, q_i, s_i, v_plot, track_cell,
    )
    t_e_spikes = np.array(t_e_spikes) if len(t_e_spikes) else np.empty(0)
    i_e_spikes = np.array(i_e_spikes, dtype=int) if len(i_e_spikes) else np.empty(0, dtype=int)
    t_i_spikes = np.array(t_i_spikes) if len(t_i_spikes) else np.empty(0)
    i_i_spikes = np.array(i_i_spikes, dtype=int) if len(i_i_spikes) else np.empty(0, dtype=int)

    return t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, v_plot


g_ee, g_ei, g_ie, g_ii = make_random_connectivity(p_ee, p_ei, p_ie, p_ii)
t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, v_plot = simulate(g_ee, g_ei, g_ie, g_ii)

if __name__ == "__main__":
    f_hat_e = round(len(t_e_spikes) / num_e / t_final * 1000)
    f_hat_i = round(len(t_i_spikes) / num_i / t_final * 1000)
    print(f"f_hat_e = {f_hat_e}")
    print(f"f_hat_i = {f_hat_i}")

    # Closeup of E-cells 1-20 (plotted at 51-70), with local firing-rate
    # maxima marked, following the matlab source's rastergram.m.

    fig, ax = plt.subplots(figsize=(8, 5))
    e_mask = (i_e_spikes < 20)
    if e_mask.any():
        ax.plot(t_e_spikes[e_mask], i_e_spikes[e_mask] + num_i + 1, '.r', markersize=6)
    ax.set_yticks([51, 70])
    ax.axis([0, t_final, 50.5, 70.5])
    ax.set_xlabel('$t$ [ms]')

    t_vec_lfp = np.arange(m_steps + 1) * dt
    rate = np.zeros(m_steps + 1)
    sigma = 3.
    # "rate" cannot be thought of as the instantaneous firing rate, but
    # as proportional to it, and since here we only want to find local
    # maxima of the instantaneous firing rate, that is good enough.
    for tk in t_e_spikes:
        rate += np.exp(-(t_vec_lfp - tk) ** 2 / (2 * sigma ** 2))
    rate_c, rate_l, rate_r = rate[1:-1], rate[:-2], rate[2:]
    ind = np.where((rate_c > rate_l) & (rate_c > rate_r))[0]
    t_peaks = ind * dt
    for tp in t_peaks:
        ax.plot([tp, tp], [50.5, 70.5], '-k', linewidth=1)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
