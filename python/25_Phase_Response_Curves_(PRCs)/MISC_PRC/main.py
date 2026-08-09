import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- WB (panel 1)


def wb_alpha_h(v):
    return 0.35 * exp(-(v + 58) / 20)


def wb_alpha_m(v):
    return 0.1 * (v + 35) / (1 - exp(-(v + 35) / 10))


def wb_alpha_n(v):
    return 0.05 * (v + 34) / (1 - exp(-0.1 * (v + 34)))


def wb_beta_h(v):
    return 5. / (exp(-0.1 * (v + 28)) + 1)


def wb_beta_m(v):
    return 4 * exp(-(v + 60) / 18)


def wb_beta_n(v):
    return 0.625 * exp(-(v + 44) / 80)


def wb_m_inf(v):
    return wb_alpha_m(v) / (wb_alpha_m(v) + wb_beta_m(v))


def wb_h_inf(v):
    return wb_alpha_h(v) / (wb_alpha_h(v) + wb_beta_h(v))


def wb_n_inf(v):
    return wb_alpha_n(v) / (wb_alpha_n(v) + wb_beta_n(v))


def wb_init(i_ext, phi_vec):
    c, g_k, g_na, g_l = 1., 9., 35., 0.1
    v_k, v_na, v_l = -90., 55., -65.
    t_final = 5000.
    dt = 0.01
    dt05 = dt / 2

    v = [-70.]
    m = [wb_m_inf(v[0])]
    h = [wb_h_inf(v[0])]
    n = [wb_n_inf(v[0])]
    t_spikes = []

    k = 0
    t = 0.
    while len(t_spikes) < 5 and t < t_final:
        vk, mk, hk, nk = v[k], m[k], h[k], n[k]
        v_inc = (g_k * nk ** 4 * (v_k - vk) + g_na * mk ** 3 * hk * (v_na - vk) + g_l * (v_l - vk) + i_ext) / c
        h_inc = wb_alpha_h(vk) * (1 - hk) - wb_beta_h(vk) * hk
        n_inc = wb_alpha_n(vk) * (1 - nk) - wb_beta_n(vk) * nk

        v_tmp = vk + dt05 * v_inc
        m_tmp = wb_m_inf(v_tmp)
        h_tmp = hk + dt05 * h_inc
        n_tmp = nk + dt05 * n_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = wb_alpha_h(v_tmp) * (1 - h_tmp) - wb_beta_h(v_tmp) * h_tmp
        n_inc = wb_alpha_n(v_tmp) * (1 - n_tmp) - wb_beta_n(v_tmp) * n_tmp

        v.append(vk + dt * v_inc)
        h.append(hk + dt * h_inc)
        n.append(nk + dt * n_inc)
        m.append(wb_m_inf(v[-1]))

        if vk >= -20 and v[-1] < -20:
            t_spike = (k * dt * (-20 - v[-1]) + (k + 1) * dt * (20 + vk)) / (vk - v[-1])
            t_spikes.append(t_spike)
        t = (k + 1) * dt
        k += 1

    num = len(phi_vec)
    T = t_spikes[4] - t_spikes[3]
    out = np.zeros((num, 3))
    for i, phi0 in enumerate(phi_vec):
        t0 = phi0 * T + t_spikes[3]
        kk = int(t0 / dt)
        frac_hi = (t0 - kk * dt) / dt
        frac_lo = ((kk + 1) * dt - t0) / dt
        out[i, 0] = v[kk + 1] * frac_hi + v[kk] * frac_lo
        out[i, 1] = h[kk + 1] * frac_hi + h[kk] * frac_lo
        out[i, 2] = n[kk + 1] * frac_hi + n[kk] * frac_lo
    return out, T


def wb_prc():
    c, g_k, g_na, g_l = 1., 9., 35., 0.1
    v_k, v_na, v_l = -90., 55., -65.
    i_ext = 0.30
    tau_r, tau_peak, tau_d = 0.5, 0.5, 2.
    tau_dq = tau_d_q_function(tau_d, tau_r, tau_peak)
    v_rev, g_syn = 0., 0.1

    N = 200
    dt = 0.01
    dt05 = dt / 2
    phi_vec = np.arange(1, N + 1) / N - 1 / (2 * N)

    initial_vector, T = wb_init(i_ext, phi_vec)
    v = initial_vector[:, 0].copy()
    m = wb_m_inf(v)
    h = initial_vector[:, 1].copy()
    n = initial_vector[:, 2].copy()
    q = np.ones(N)
    s = np.zeros(N)

    t_star = np.full(N, np.nan)
    num_spikes = np.zeros(N, dtype=int)

    k = 0
    while np.min(num_spikes) < 1:
        k += 1

        v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v)
                  + g_syn * s * (v_rev - v) + i_ext) / c
        h_inc = wb_alpha_h(v) * (1 - h) - wb_beta_h(v) * h
        n_inc = wb_alpha_n(v) * (1 - n) - wb_beta_n(v) * n
        q_inc = -q / tau_dq
        s_inc = q * (1 - s) / tau_r - s / tau_d

        v_tmp = v + dt05 * v_inc
        m_tmp = wb_m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        q_tmp = q + dt05 * q_inc
        s_tmp = s + dt05 * s_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + g_syn * s_tmp * (v_rev - v_tmp) + i_ext) / c
        h_inc = wb_alpha_h(v_tmp) * (1 - h_tmp) - wb_beta_h(v_tmp) * h_tmp
        n_inc = wb_alpha_n(v_tmp) * (1 - n_tmp) - wb_beta_n(v_tmp) * n_tmp
        q_inc = -q_tmp / tau_dq
        s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d

        v_old = v.copy()
        v = v + dt * v_inc
        m = wb_m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc
        q = q + dt * q_inc
        s = s + dt * s_inc

        ind = np.where((v_old >= -20) & (v < -20))[0]
        for i in ind:
            if num_spikes[i] == 0:
                t_star[i] = ((k - 1) * dt * (-20 - v[i]) + k * dt * (v_old[i] + 20)) / (v_old[i] - v[i])
            num_spikes[i] += 1

    g_vec = -t_star / T + 1 - phi_vec
    return phi_vec, g_vec


# ---------------------------------------------------------------- HH (panel 2)


def hh_alpha_h(v):
    return 0.07 * exp(-(v + 70) / 20)


def hh_alpha_m(v):
    return np.where(np.abs(v + 45) > 1e-8, (v + 45) / 10 / (1 - exp(-(v + 45) / 10)), 1.)


def hh_alpha_n(v):
    return 0.01 * (-60. - v) / (exp((-60 - v) / 10) - 1)


def hh_beta_h(v):
    return 1. / (exp(-(v + 40) / 10) + 1)


def hh_beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def hh_beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def hh_m_inf(v):
    return hh_alpha_m(v) / (hh_alpha_m(v) + hh_beta_m(v))


def hh_h_inf(v):
    return hh_alpha_h(v) / (hh_alpha_h(v) + hh_beta_h(v))


def hh_n_inf(v):
    return hh_alpha_n(v) / (hh_alpha_n(v) + hh_beta_n(v))


def hh_init(i_ext, phi_vec):
    c, g_k, g_na, g_l = 1., 36., 120., 0.3
    v_k, v_na, v_l = -82., 45., -59.
    t_final = 5000.
    dt = 0.005
    dt05 = dt / 2

    v = [-70.]
    m = [hh_m_inf(v[0])]
    h = [hh_h_inf(v[0])]
    n = [hh_n_inf(v[0])]
    t_spikes = []

    k = 0
    t = 0.
    while len(t_spikes) < 5 and t < t_final:
        vk, mk, hk, nk = v[k], m[k], h[k], n[k]
        v_inc = (g_k * nk ** 4 * (v_k - vk) + g_na * mk ** 3 * hk * (v_na - vk) + g_l * (v_l - vk) + i_ext) / c
        m_inc = hh_alpha_m(vk) * (1 - mk) - hh_beta_m(vk) * mk
        h_inc = hh_alpha_h(vk) * (1 - hk) - hh_beta_h(vk) * hk
        n_inc = hh_alpha_n(vk) * (1 - nk) - hh_beta_n(vk) * nk

        v_tmp = vk + dt05 * v_inc
        m_tmp = mk + dt05 * m_inc
        h_tmp = hk + dt05 * h_inc
        n_tmp = nk + dt05 * n_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        m_inc = hh_alpha_m(v_tmp) * (1 - m_tmp) - hh_beta_m(v_tmp) * m_tmp
        h_inc = hh_alpha_h(v_tmp) * (1 - h_tmp) - hh_beta_h(v_tmp) * h_tmp
        n_inc = hh_alpha_n(v_tmp) * (1 - n_tmp) - hh_beta_n(v_tmp) * n_tmp

        v.append(vk + dt * v_inc)
        m.append(mk + dt * m_inc)
        h.append(hk + dt * h_inc)
        n.append(nk + dt * n_inc)

        if vk >= -20 and v[-1] < -20:
            t_spike = (k * dt * (-20 - v[-1]) + (k + 1) * dt * (20 + vk)) / (vk - v[-1])
            t_spikes.append(t_spike)
        t = (k + 1) * dt
        k += 1

    num = len(phi_vec)
    T = t_spikes[4] - t_spikes[3]
    out = np.zeros((num, 4))
    for i, phi0 in enumerate(phi_vec):
        t0 = phi0 * T + t_spikes[3]
        kk = int(t0 / dt)
        frac_hi = (t0 - kk * dt) / dt
        frac_lo = ((kk + 1) * dt - t0) / dt
        out[i, 0] = v[kk + 1] * frac_hi + v[kk] * frac_lo
        out[i, 1] = m[kk + 1] * frac_hi + m[kk] * frac_lo
        out[i, 2] = h[kk + 1] * frac_hi + h[kk] * frac_lo
        out[i, 3] = n[kk + 1] * frac_hi + n[kk] * frac_lo
    return out, T


def hh_prc():
    c, g_k, g_na, g_l = 1., 36., 120., 0.3
    v_k, v_na, v_l = -82., 45., -59.
    i_ext = 10.
    tau_r, tau_peak, tau_d = 0.5, 0.5, 2.
    tau_dq = tau_d_q_function(tau_d, tau_r, tau_peak)
    v_rev, g_syn = 0., 0.1

    N = 500
    dt = 0.005
    dt05 = dt / 2
    phi_vec = np.arange(1, N + 1) / N - 1 / (2 * N)

    initial_vector, T = hh_init(i_ext, phi_vec)
    v = initial_vector[:, 0].copy()
    m = initial_vector[:, 1].copy()
    h = initial_vector[:, 2].copy()
    n = initial_vector[:, 3].copy()
    q = np.ones(N)
    s = np.zeros(N)

    t_star = np.full(N, np.nan)
    num_spikes = np.zeros(N, dtype=int)

    k = 0
    while np.min(num_spikes) < 1:
        k += 1

        v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v)
                  + g_syn * s * (v_rev - v) + i_ext) / c
        m_inc = hh_alpha_m(v) * (1 - m) - hh_beta_m(v) * m
        h_inc = hh_alpha_h(v) * (1 - h) - hh_beta_h(v) * h
        n_inc = hh_alpha_n(v) * (1 - n) - hh_beta_n(v) * n
        q_inc = -q / tau_dq
        s_inc = q * (1 - s) / tau_r - s / tau_d

        v_tmp = v + dt05 * v_inc
        m_tmp = m + dt05 * m_inc
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        q_tmp = q + dt05 * q_inc
        s_tmp = s + dt05 * s_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + g_syn * s_tmp * (v_rev - v_tmp) + i_ext) / c
        m_inc = hh_alpha_m(v_tmp) * (1 - m_tmp) - hh_beta_m(v_tmp) * m_tmp
        h_inc = hh_alpha_h(v_tmp) * (1 - h_tmp) - hh_beta_h(v_tmp) * h_tmp
        n_inc = hh_alpha_n(v_tmp) * (1 - n_tmp) - hh_beta_n(v_tmp) * n_tmp
        q_inc = -q_tmp / tau_dq
        s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d

        v_old = v.copy()
        v = v + dt * v_inc
        m = m + dt * m_inc
        h = h + dt * h_inc
        n = n + dt * n_inc
        q = q + dt * q_inc
        s = s + dt * s_inc

        ind = np.where((v_old >= -20) & (v < -20))[0]
        for i in ind:
            if num_spikes[i] == 0:
                t_star[i] = ((k - 1) * dt * (-20 - v[i]) + k * dt * (v_old[i] + 20)) / (v_old[i] - v[i])
            num_spikes[i] += 1

    g_vec = -t_star / T + 1 - phi_vec
    return phi_vec, g_vec


# ------------------------------------------------------------ Erisir (panel 3)


def er_alpha_h(v):
    return 0.0035 / exp(v / 24.186)


def er_alpha_m(v):
    return 40 * (75.5 - v) / (exp((75.5 - v) / 13.5) - 1)


def er_alpha_n(v):
    return (95 - v) / (exp((95 - v) / 11.8) - 1)


def er_beta_h(v):
    return -0.017 * (v + 51.25) / (exp(-(v + 51.25) / 5.2) - 1)


def er_beta_m(v):
    return 1.2262 / exp(v / 42.248)


def er_beta_n(v):
    return 0.025 / exp(v / 22.222)


def er_m_inf(v):
    return er_alpha_m(v) / (er_alpha_m(v) + er_beta_m(v))


def er_h_inf(v):
    return er_alpha_h(v) / (er_alpha_h(v) + er_beta_h(v))


def er_n_inf(v):
    return er_alpha_n(v) / (er_alpha_n(v) + er_beta_n(v))


def erisir_init(i_ext, phi_vec):
    c, g_k, g_na, g_l = 1., 224., 112., 0.5
    v_k, v_na, v_l = -90., 60., -70.
    t_final = 5000.
    dt = 0.002
    dt05 = dt / 2

    v = [-70.]
    m = [er_m_inf(v[0])]
    h = [er_h_inf(v[0])]
    n = [er_n_inf(v[0])]
    t_spikes = []

    k = 0
    t = 0.
    while len(t_spikes) < 5 and t < t_final:
        vk, mk, hk, nk = v[k], m[k], h[k], n[k]
        v_inc = (g_k * nk ** 2 * (v_k - vk) + g_na * mk ** 3 * hk * (v_na - vk) + g_l * (v_l - vk) + i_ext) / c
        h_inc = er_alpha_h(vk) * (1 - hk) - er_beta_h(vk) * hk
        n_inc = er_alpha_n(vk) * (1 - nk) - er_beta_n(vk) * nk

        v_tmp = vk + dt05 * v_inc
        m_tmp = er_m_inf(v_tmp)
        h_tmp = hk + dt05 * h_inc
        n_tmp = nk + dt05 * n_inc

        v_inc = (g_k * n_tmp ** 2 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = er_alpha_h(v_tmp) * (1 - h_tmp) - er_beta_h(v_tmp) * h_tmp
        n_inc = er_alpha_n(v_tmp) * (1 - n_tmp) - er_beta_n(v_tmp) * n_tmp

        v.append(vk + dt * v_inc)
        h.append(hk + dt * h_inc)
        n.append(nk + dt * n_inc)
        m.append(er_m_inf(v[-1]))

        if vk >= -20 and v[-1] < -20:
            t_spike = (k * dt * (-20 - v[-1]) + (k + 1) * dt * (20 + vk)) / (vk - v[-1])
            t_spikes.append(t_spike)
        t = (k + 1) * dt
        k += 1

    num = len(phi_vec)
    T = t_spikes[4] - t_spikes[3]
    out = np.zeros((num, 3))
    for i, phi0 in enumerate(phi_vec):
        t0 = phi0 * T + t_spikes[3]
        kk = int(t0 / dt)
        frac_hi = (t0 - kk * dt) / dt
        frac_lo = ((kk + 1) * dt - t0) / dt
        out[i, 0] = v[kk + 1] * frac_hi + v[kk] * frac_lo
        out[i, 1] = h[kk + 1] * frac_hi + h[kk] * frac_lo
        out[i, 2] = n[kk + 1] * frac_hi + n[kk] * frac_lo
    return out, T


def erisir_prc():
    c, g_k, g_na, g_l = 1., 224., 112., 0.5
    v_k, v_na, v_l = -90., 60., -70.
    i_ext = 7.1
    tau_r, tau_peak, tau_d = 0.5, 0.5, 2.
    tau_dq = tau_d_q_function(tau_d, tau_r, tau_peak)
    v_rev, g_syn = 0., 0.1

    N = 200
    dt = 0.002
    dt05 = dt / 2
    phi_vec = np.arange(1, N + 1) / N - 1 / (2 * N)

    initial_vector, T = erisir_init(i_ext, phi_vec)
    v = initial_vector[:, 0].copy()
    m = er_m_inf(v)
    h = initial_vector[:, 1].copy()
    n = initial_vector[:, 2].copy()
    q = np.ones(N)
    s = np.zeros(N)

    t_star = np.full(N, np.nan)
    num_spikes = np.zeros(N, dtype=int)

    k = 0
    while np.min(num_spikes) < 1:
        k += 1

        v_inc = (g_k * n ** 2 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v)
                  + g_syn * s * (v_rev - v) + i_ext) / c
        h_inc = er_alpha_h(v) * (1 - h) - er_beta_h(v) * h
        n_inc = er_alpha_n(v) * (1 - n) - er_beta_n(v) * n
        q_inc = -q / tau_dq
        s_inc = q * (1 - s) / tau_r - s / tau_d

        v_tmp = v + dt05 * v_inc
        m_tmp = er_m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        q_tmp = q + dt05 * q_inc
        s_tmp = s + dt05 * s_inc

        v_inc = (g_k * n_tmp ** 2 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + g_syn * s_tmp * (v_rev - v_tmp) + i_ext) / c
        h_inc = er_alpha_h(v_tmp) * (1 - h_tmp) - er_beta_h(v_tmp) * h_tmp
        n_inc = er_alpha_n(v_tmp) * (1 - n_tmp) - er_beta_n(v_tmp) * n_tmp
        q_inc = -q_tmp / tau_dq
        s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d

        v_old = v.copy()
        v = v + dt * v_inc
        m = er_m_inf(v)
        h = h + dt * h_inc
        n = n + dt * n_inc
        q = q + dt * q_inc
        s = s + dt * s_inc

        ind = np.where((v_old >= -20) & (v < -20))[0]
        for i in ind:
            if num_spikes[i] == 0:
                t_star[i] = ((k - 1) * dt * (-20 - v[i]) + k * dt * (v_old[i] + 20)) / (v_old[i] - v[i])
            num_spikes[i] += 1

    g_vec = -t_star / T + 1 - phi_vec
    return phi_vec, g_vec


# ------------------------------------------------------------------- shared


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


phi_vec_wb, g_vec_wb = wb_prc()
phi_vec_hh, g_vec_hh = hh_prc()
phi_vec_erisir, g_vec_erisir = erisir_prc()

if __name__ == "__main__":

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

    axes[0].plot(phi_vec_wb, g_vec_wb, '-k', linewidth=2)
    axes[0].plot([0, 1], [1, 0], '--k', linewidth=1)
    axes[0].axis([0, 1, 0, 1])
    axes[0].set_box_aspect(1)
    axes[0].set_xlabel(r'$\varphi$')
    axes[0].set_ylabel('$g$')
    axes[0].set_title('WB')

    axes[1].plot(phi_vec_hh, g_vec_hh, '-k', linewidth=2)
    axes[1].axis([0, 1, -0.1, 0.1])
    axes[1].set_box_aspect(1)
    axes[1].set_xlabel(r'$\varphi$')
    axes[1].set_title('HH')

    axes[2].plot(phi_vec_erisir, g_vec_erisir, '-k', linewidth=2)
    axes[2].axis([0, 1, 0, 0.19])
    axes[2].set_box_aspect(1)
    axes[2].set_xlabel(r'$\varphi$')
    axes[2].set_title('Erisir')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
