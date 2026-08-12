"""Neuron and synapse equations shared by the Chapter 38 examples."""

import math

import numpy as np
from numba import njit
from numba.typed import List


def _ratio(x, scale):
    """Return x / (1 - exp(-x/scale)), including its value at zero."""
    x = np.asarray(x)
    return np.divide(x, -np.expm1(-x / scale), out=np.full_like(x, scale), where=x != 0)


def m_e_inf(v):
    v = np.asarray(v)
    alpha = 0.32 * _ratio(v + 54, 4)
    beta = 0.28 * np.divide(
        v + 27,
        np.expm1((v + 27) / 5),
        out=np.full_like(v, 5.0),
        where=v != -27,
    )
    return alpha / (alpha + beta)


def h_e_inf(v):
    alpha = 0.128 * np.exp(-(np.asarray(v) + 50) / 18)
    beta = 4 / (1 + np.exp(-(np.asarray(v) + 27) / 5))
    return alpha / (alpha + beta)


def tau_h_e(v):
    alpha = 0.128 * np.exp(-(np.asarray(v) + 50) / 18)
    beta = 4 / (1 + np.exp(-(np.asarray(v) + 27) / 5))
    return 1 / (alpha + beta)


def n_e_inf(v):
    v = np.asarray(v)
    alpha = 0.032 * _ratio(v + 52, 5)
    beta = 0.5 * np.exp(-(v + 57) / 40)
    return alpha / (alpha + beta)


def tau_n_e(v):
    v = np.asarray(v)
    alpha = 0.032 * _ratio(v + 52, 5)
    beta = 0.5 * np.exp(-(v + 57) / 40)
    return 1 / (alpha + beta)


def m_i_inf(v):
    v = np.asarray(v)
    alpha = 0.1 * _ratio(v + 35, 10)
    beta = 4 * np.exp(-(v + 60) / 18)
    return alpha / (alpha + beta)


def h_i_inf(v):
    alpha = 0.07 * np.exp(-(np.asarray(v) + 58) / 20)
    beta = 1 / (np.exp(-0.1 * (np.asarray(v) + 28)) + 1)
    return alpha / (alpha + beta)


def tau_h_i(v):
    alpha = 0.07 * np.exp(-(np.asarray(v) + 58) / 20)
    beta = 1 / (np.exp(-0.1 * (np.asarray(v) + 28)) + 1)
    return 1 / (alpha + beta) / 5


def n_i_inf(v):
    v = np.asarray(v)
    alpha = 0.01 * _ratio(v + 34, 10)
    beta = 0.125 * np.exp(-(v + 44) / 80)
    return alpha / (alpha + beta)


def tau_n_i(v):
    v = np.asarray(v)
    alpha = 0.01 * _ratio(v + 34, 10)
    beta = 0.125 * np.exp(-(v + 44) / 80)
    return 1 / (alpha + beta) / 5


# ---------------------------------------------------------------------------
# Scalar njit mirrors of the gating functions above, used only by the
# _poisson_ping_loop hot loop below (numba needs scalar math.exp/expm1
# calls, not vectorized numpy, to compile a tight per-neuron loop).


@njit
def _ratio_s(x, scale):
    if x == 0.0:
        return scale
    return x / (-math.expm1(-x / scale))


@njit
def _m_e_inf_s(v):
    alpha = 0.32 * _ratio_s(v + 54, 4.0)
    if v + 27 == 0.0:
        beta = 0.28 * 5.0
    else:
        beta = 0.28 * (v + 27) / math.expm1((v + 27) / 5)
    return alpha / (alpha + beta)


@njit
def _h_e_inf_s(v):
    alpha = 0.128 * math.exp(-(v + 50) / 18)
    beta = 4.0 / (1 + math.exp(-(v + 27) / 5))
    return alpha / (alpha + beta)


@njit
def _tau_h_e_s(v):
    alpha = 0.128 * math.exp(-(v + 50) / 18)
    beta = 4.0 / (1 + math.exp(-(v + 27) / 5))
    return 1.0 / (alpha + beta)


@njit
def _n_e_inf_s(v):
    alpha = 0.032 * _ratio_s(v + 52, 5.0)
    beta = 0.5 * math.exp(-(v + 57) / 40)
    return alpha / (alpha + beta)


@njit
def _tau_n_e_s(v):
    alpha = 0.032 * _ratio_s(v + 52, 5.0)
    beta = 0.5 * math.exp(-(v + 57) / 40)
    return 1.0 / (alpha + beta)


@njit
def _m_i_inf_s(v):
    alpha = 0.1 * _ratio_s(v + 35, 10.0)
    beta = 4.0 * math.exp(-(v + 60) / 18)
    return alpha / (alpha + beta)


@njit
def _h_i_inf_s(v):
    alpha = 0.07 * math.exp(-(v + 58) / 20)
    beta = 1.0 / (math.exp(-0.1 * (v + 28)) + 1)
    return alpha / (alpha + beta)


@njit
def _tau_h_i_s(v):
    alpha = 0.07 * math.exp(-(v + 58) / 20)
    beta = 1.0 / (math.exp(-0.1 * (v + 28)) + 1)
    return 1.0 / (alpha + beta) / 5.0


@njit
def _n_i_inf_s(v):
    alpha = 0.01 * _ratio_s(v + 34, 10.0)
    beta = 0.125 * math.exp(-(v + 44) / 80)
    return alpha / (alpha + beta)


@njit
def _tau_n_i_s(v):
    alpha = 0.01 * _ratio_s(v + 34, 10.0)
    beta = 0.125 * math.exp(-(v + 44) / 80)
    return 1.0 / (alpha + beta) / 5.0


@njit
def _poisson_ping_loop(steps, dt, num_e, num_i, tdqe, tdqi, mu,
                        has_pulse, pulse_period, pulse_phase,
                        ve, he, ne, qe, se, q_stoch, s_stoch,
                        vi, hi, ni, qi, si, poisson_draws):
    e_times = List.empty_list(np.float64)
    e_indices = List.empty_list(np.int64)
    i_times = List.empty_list(np.float64)
    i_indices = List.empty_list(np.int64)

    dve = np.empty(num_e); dhe = np.empty(num_e); dne = np.empty(num_e)
    dqe = np.empty(num_e); dse = np.empty(num_e); dqs = np.empty(num_e); dss = np.empty(num_e)
    dvi = np.empty(num_i); dhi = np.empty(num_i); dni = np.empty(num_i)
    dqi = np.empty(num_i); dsi = np.empty(num_i)

    ve_m = np.empty(num_e); he_m = np.empty(num_e); ne_m = np.empty(num_e)
    qe_m = np.empty(num_e); se_m = np.empty(num_e); qs_m = np.empty(num_e); ss_m = np.empty(num_e)
    vi_m = np.empty(num_i); hi_m = np.empty(num_i); ni_m = np.empty(num_i)
    qi_m = np.empty(num_i); si_m = np.empty(num_i)

    old_ve = np.empty(num_e)
    old_vi = np.empty(num_i)

    for k in range(steps):
        for j in range(num_e):
            old_ve[j] = ve[j]
        for j in range(num_i):
            old_vi[j] = vi[j]

        if has_pulse:
            phase = k * dt / pulse_period - pulse_phase
            pulse_term = 0.1 * (math.exp(4 * math.cos(math.pi * phase) ** 2) - 1) / mu
        else:
            pulse_term = 0.0

        se_mean = 0.0
        for j in range(num_e):
            se_mean += se[j]
        se_mean /= num_e
        si_mean = 0.0
        for j in range(num_i):
            si_mean += si[j]
        si_mean /= num_i

        for j in range(num_e):
            v = ve[j]
            drive = 0.6 + 0.1 * s_stoch[j] * (0.0 - v)
            if has_pulse and j < 5:
                drive += pulse_term
            me = _m_e_inf_s(v)
            dve[j] = (0.1 * (-67 - v) + 80 * ne[j] ** 4 * (-100 - v)
                      + 100 * me ** 3 * he[j] * (50 - v)
                      + 1.25 * si_mean * (-75 - v) + drive)
            dhe[j] = (_h_e_inf_s(v) - he[j]) / _tau_h_e_s(v)
            dne[j] = (_n_e_inf_s(v) - ne[j]) / _tau_n_e_s(v)
            th = math.tanh(v / 10)
            dqe[j] = (1 + th) / 2 * (1 - qe[j]) / 0.1 - qe[j] / tdqe
            dse[j] = qe[j] * (1 - se[j]) / 0.3 - se[j] / 3
            dqs[j] = -q_stoch[j] / tdqe
            dss[j] = q_stoch[j] * (1 - s_stoch[j]) / 0.3 - s_stoch[j] / 3

        for j in range(num_i):
            v = vi[j]
            mi = _m_i_inf_s(v)
            dvi[j] = (0.1 * (-65 - v) + 9 * ni[j] ** 4 * (-90 - v)
                      + 35 * mi ** 3 * hi[j] * (55 - v)
                      + 1.25 * se_mean * (0.0 - v) + 0.4 * si_mean * (-75 - v) + 0.6)
            dhi[j] = (_h_i_inf_s(v) - hi[j]) / _tau_h_i_s(v)
            dni[j] = (_n_i_inf_s(v) - ni[j]) / _tau_n_i_s(v)
            thi = math.tanh(v / 10)
            dqi[j] = (1 + thi) / 2 * (1 - qi[j]) / 0.1 - qi[j] / tdqi
            dsi[j] = qi[j] * (1 - si[j]) / 0.3 - si[j] / 9

        for j in range(num_e):
            ve_m[j] = ve[j] + 0.5 * dt * dve[j]
            he_m[j] = he[j] + 0.5 * dt * dhe[j]
            ne_m[j] = ne[j] + 0.5 * dt * dne[j]
            qe_m[j] = qe[j] + 0.5 * dt * dqe[j]
            se_m[j] = se[j] + 0.5 * dt * dse[j]
            qs_m[j] = q_stoch[j] + 0.5 * dt * dqs[j]
            ss_m[j] = s_stoch[j] + 0.5 * dt * dss[j]
        for j in range(num_i):
            vi_m[j] = vi[j] + 0.5 * dt * dvi[j]
            hi_m[j] = hi[j] + 0.5 * dt * dhi[j]
            ni_m[j] = ni[j] + 0.5 * dt * dni[j]
            qi_m[j] = qi[j] + 0.5 * dt * dqi[j]
            si_m[j] = si[j] + 0.5 * dt * dsi[j]

        if has_pulse:
            phase_mid = (k + 0.5) * dt / pulse_period - pulse_phase
            pulse_term_mid = 0.1 * (math.exp(4 * math.cos(math.pi * phase_mid) ** 2) - 1) / mu
        else:
            pulse_term_mid = 0.0

        se_mean_m = 0.0
        for j in range(num_e):
            se_mean_m += se_m[j]
        se_mean_m /= num_e
        si_mean_m = 0.0
        for j in range(num_i):
            si_mean_m += si_m[j]
        si_mean_m /= num_i

        for j in range(num_e):
            v = ve_m[j]
            drive = 0.6 + 0.1 * ss_m[j] * (0.0 - v)
            if has_pulse and j < 5:
                drive += pulse_term_mid
            me = _m_e_inf_s(v)
            dve[j] = (0.1 * (-67 - v) + 80 * ne_m[j] ** 4 * (-100 - v)
                      + 100 * me ** 3 * he_m[j] * (50 - v)
                      + 1.25 * si_mean_m * (-75 - v) + drive)
            dhe[j] = (_h_e_inf_s(v) - he_m[j]) / _tau_h_e_s(v)
            dne[j] = (_n_e_inf_s(v) - ne_m[j]) / _tau_n_e_s(v)
            th = math.tanh(v / 10)
            dqe[j] = (1 + th) / 2 * (1 - qe_m[j]) / 0.1 - qe_m[j] / tdqe
            dse[j] = qe_m[j] * (1 - se_m[j]) / 0.3 - se_m[j] / 3
            dqs[j] = -qs_m[j] / tdqe
            dss[j] = qs_m[j] * (1 - ss_m[j]) / 0.3 - ss_m[j] / 3

        for j in range(num_i):
            v = vi_m[j]
            mi = _m_i_inf_s(v)
            dvi[j] = (0.1 * (-65 - v) + 9 * ni_m[j] ** 4 * (-90 - v)
                      + 35 * mi ** 3 * hi_m[j] * (55 - v)
                      + 1.25 * se_mean_m * (0.0 - v) + 0.4 * si_mean_m * (-75 - v) + 0.6)
            dhi[j] = (_h_i_inf_s(v) - hi_m[j]) / _tau_h_i_s(v)
            dni[j] = (_n_i_inf_s(v) - ni_m[j]) / _tau_n_i_s(v)
            thi = math.tanh(v / 10)
            dqi[j] = (1 + thi) / 2 * (1 - qi_m[j]) / 0.1 - qi_m[j] / tdqi
            dsi[j] = qi_m[j] * (1 - si_m[j]) / 0.3 - si_m[j] / 9

        for j in range(num_e):
            ve[j] = ve[j] + dt * dve[j]
            he[j] = he[j] + dt * dhe[j]
            ne[j] = ne[j] + dt * dne[j]
            qe[j] = qe[j] + dt * dqe[j]
            se[j] = se[j] + dt * dse[j]
            q_stoch[j] = q_stoch[j] + dt * dqs[j]
            s_stoch[j] = s_stoch[j] + dt * dss[j]
        for j in range(num_i):
            vi[j] = vi[j] + dt * dvi[j]
            hi[j] = hi[j] + dt * dhi[j]
            ni[j] = ni[j] + dt * dni[j]
            qi[j] = qi[j] + dt * dqi[j]
            si[j] = si[j] + dt * dsi[j]

        for j in range(num_e):
            if old_ve[j] > -20 and ve[j] <= -20:
                e_indices.append(j)
                e_times.append(((-20 - ve[j]) * k * dt + (old_ve[j] + 20) * (k + 1) * dt)
                                / (old_ve[j] - ve[j]))
        for j in range(num_i):
            if old_vi[j] > -20 and vi[j] <= -20:
                i_indices.append(j)
                i_times.append(((-20 - vi[j]) * k * dt + (old_vi[j] + 20) * (k + 1) * dt)
                                / (old_vi[j] - vi[j]))

        for j in range(num_e):
            if poisson_draws[k, j]:
                q_stoch[j] = 1.0

    return e_times, e_indices, i_times, i_indices


def _tau_peak(tau_d, tau_r, tau_dq):
    dt = 0.01
    s = t = 0.0
    derivative = 1 / tau_r
    while derivative > 0:
        t_old, derivative_old = t, derivative
        s_mid = s + 0.5 * dt * derivative
        derivative_mid = np.exp(-(t + 0.5 * dt) / tau_dq) * (1 - s_mid) / tau_r - s_mid / tau_d
        s += dt * derivative_mid
        t += dt
        derivative = np.exp(-t / tau_dq) * (1 - s) / tau_r - s / tau_d
    return (t_old * -derivative + t * derivative_old) / (derivative_old - derivative)


def tau_d_q(tau_d, tau_r, tau_peak):
    left = 1.0
    while _tau_peak(tau_d, tau_r, left) > tau_peak:
        left /= 2
    right = tau_r
    while _tau_peak(tau_d, tau_r, right) < tau_peak:
        right *= 2
    while right - left > 1e-12:
        middle = (left + right) / 2
        if _tau_peak(tau_d, tau_r, middle) <= tau_peak:
            left = middle
        else:
            right = middle
    return (left + right) / 2


def simulate_two_cell(g_ei, t_final=200.0, mean_inhibition=None):
    """Simulate the driven RTM/WB pair used by GAMMA_COHERENCE_1 and 2."""
    dt = 0.01
    steps = round(t_final / dt)
    t = np.arange(steps + 1) * dt
    sample = np.arange(1, 3001) / 3000
    main = 0.5 * (np.exp(20 * np.cos(np.pi * t / 25) ** 2) - 1)
    main /= np.mean(np.exp(20 * np.cos(np.pi * sample) ** 2) - 1)
    distractor = 0.5 * (np.exp(0.1 * np.cos(np.pi * t / 40) ** 2) - 1)
    distractor /= np.mean(np.exp(0.1 * np.cos(np.pi * sample) ** 2) - 1)

    tdqe, tdqi = tau_d_q(3, 0.5, 0.5), tau_d_q(9, 0.5, 0.5)
    ve, vi = -70.0, -75.0
    he, ne = float(h_e_inf(ve)), float(n_e_inf(ve))
    hi, ni = float(h_i_inf(vi)), float(n_i_inf(vi))
    qe = se = qi = si = 0.0
    ve_trace, vi_trace, si_trace = np.empty(steps + 1), np.empty(steps + 1), np.empty(steps + 1)
    ve_trace[0], vi_trace[0], si_trace[0] = ve, vi, si
    e_spikes, i_spikes = [], []

    def slopes(ve_, he_, ne_, qe_, se_, vi_, hi_, ni_, qi_, si_, drive, inhibition):
        me, mi = m_e_inf(ve_), m_i_inf(vi_)
        dve = (0.1 * (-67 - ve_) + 80 * ne_ ** 4 * (-100 - ve_)
               + 100 * me ** 3 * he_ * (50 - ve_) + 0.5 * inhibition * (-75 - ve_) + drive)
        dvi = (0.1 * (-65 - vi_) + 9 * ni_ ** 4 * (-90 - vi_)
               + 35 * mi ** 3 * hi_ * (55 - vi_) + g_ei * se_ * (0 - vi_)
               + 0.15 * si_ * (-75 - vi_))
        return np.array([
            dve, (h_e_inf(ve_) - he_) / tau_h_e(ve_), (n_e_inf(ve_) - ne_) / tau_n_e(ve_),
            (1 + np.tanh(ve_ / 10)) / 2 * (1 - qe_) / 0.1 - qe_ / tdqe,
            qe_ * (1 - se_) / 0.5 - se_ / 3,
            dvi, (h_i_inf(vi_) - hi_) / tau_h_i(vi_), (n_i_inf(vi_) - ni_) / tau_n_i(vi_),
            (1 + np.tanh(vi_ / 10)) / 2 * (1 - qi_) / 0.1 - qi_ / tdqi,
            qi_ * (1 - si_) / 0.5 - si_ / 9,
        ], dtype=float)

    state = np.array([ve, he, ne, qe, se, vi, hi, ni, qi, si])
    for k in range(steps):
        ve_old, vi_old = state[0], state[5]
        inhibition = state[9] if mean_inhibition is None else mean_inhibition
        first = slopes(*state, main[k] + distractor[k], inhibition)
        middle = state + 0.5 * dt * first
        inhibition_mid = middle[9] if mean_inhibition is None else mean_inhibition
        second = slopes(*middle, (main[k] + main[k + 1] + distractor[k] + distractor[k + 1]) / 2, inhibition_mid)
        state += dt * second
        ve, vi = state[0], state[5]
        if ve_old > -20 >= ve:
            e_spikes.append((t[k] * (-20 - ve) + t[k + 1] * (ve_old + 20)) / (ve_old - ve))
        if vi_old > -20 >= vi:
            i_spikes.append((t[k] * (-20 - vi) + t[k + 1] * (vi_old + 20)) / (vi_old - vi))
        ve_trace[k + 1], vi_trace[k + 1], si_trace[k + 1] = ve, vi, state[9]

    return {
        "t": t, "v_e": ve_trace, "v_i": vi_trace, "s_i": si_trace,
        "t_e_spikes": np.asarray(e_spikes), "t_i_spikes": np.asarray(i_spikes),
        "i_main": main, "i_dist": distractor,
    }


def simulate_poisson_ping(seed=63806, t_final=500.0, pulse_period=None, pulse_phase=0.0):
    """Simulate the Chapter 38 all-to-all PING network."""
    rng = np.random.default_rng(seed)
    num_e, num_i, dt = 200, 50, 0.01
    steps = round(t_final / dt)
    tdqe, tdqi = tau_d_q(3, 0.3, 0.3), tau_d_q(9, 0.3, 0.3)

    # Random points near each population's repetitive-firing orbit provide
    # the randomized phases used by the MATLAB rtm_init/wb_init helpers.
    ve = rng.uniform(-70, -45, num_e)
    vi = rng.uniform(-70, -45, num_i)
    he, ne, hi, ni = h_e_inf(ve), n_e_inf(ve), h_i_inf(vi), n_i_inf(vi)
    qe = np.zeros(num_e); se = np.zeros(num_e)
    qi = np.zeros(num_i); si = np.zeros(num_i)
    q_stoch = np.zeros(num_e); s_stoch = np.zeros(num_e)

    mu = np.mean(np.exp(4 * np.cos(np.pi * np.arange(1, 1001) / 100) ** 2) - 1)

    # rng.random((steps, num_e)) draws from the same underlying bit stream,
    # in the same order, as calling rng.random(num_e) once per step -- this
    # lets the per-step Poisson trigger be precomputed outside the hot loop
    # while keeping the exact reproducibility guaranteed by the seed.
    poisson_draws = rng.random((steps, num_e)) < 40 / 1000 * dt

    has_pulse = pulse_period is not None
    e_times, e_indices, i_times, i_indices = _poisson_ping_loop(
        steps, dt, num_e, num_i, tdqe, tdqi, mu,
        has_pulse, pulse_period if has_pulse else 1.0, pulse_phase,
        ve, he, ne, qe, se, q_stoch, s_stoch, vi, hi, ni, qi, si, poisson_draws,
    )
    e_times = np.array(e_times) if len(e_times) else np.empty(0)
    e_indices = np.array(e_indices, dtype=int) if len(e_indices) else np.empty(0, dtype=int)
    i_times = np.array(i_times) if len(i_times) else np.empty(0)
    i_indices = np.array(i_indices, dtype=int) if len(i_indices) else np.empty(0, dtype=int)

    return {
        "t_e_spikes": e_times, "i_e_spikes": e_indices,
        "t_i_spikes": i_times, "i_i_spikes": i_indices,
        "f_hat_e": len(e_times) / num_e / t_final * 1000,
        "f_hat_i": len(i_times) / num_i / t_final * 1000,
    }
