"""Neuron and synapse equations shared by the Chapter 38 examples."""

import numpy as np


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
    e_times, e_indices, i_times, i_indices = [], [], [], []

    mu = np.mean(np.exp(4 * np.cos(np.pi * np.arange(1, 1001) / 100) ** 2) - 1)

    def slopes(ve_, he_, ne_, qe_, se_, qs_, ss_, vi_, hi_, ni_, qi_, si_, pulse):
        drive_e = 0.6 + 0.1 * ss_ * (0 - ve_)
        if pulse_period is not None:
            drive_e = drive_e.copy()
            drive_e[:5] += 0.1 * (np.exp(4 * np.cos(np.pi * pulse) ** 2) - 1) / mu
        dve = (0.1 * (-67 - ve_) + 80 * ne_ ** 4 * (-100 - ve_)
               + 100 * m_e_inf(ve_) ** 3 * he_ * (50 - ve_)
               + 1.25 * si_.mean() * (-75 - ve_) + drive_e)
        dvi = (0.1 * (-65 - vi_) + 9 * ni_ ** 4 * (-90 - vi_)
               + 35 * m_i_inf(vi_) ** 3 * hi_ * (55 - vi_)
               + 1.25 * se_.mean() * (0 - vi_) + 0.4 * si_.mean() * (-75 - vi_) + 0.6)
        return (
            dve, (h_e_inf(ve_) - he_) / tau_h_e(ve_), (n_e_inf(ve_) - ne_) / tau_n_e(ve_),
            (1 + np.tanh(ve_ / 10)) / 2 * (1 - qe_) / 0.1 - qe_ / tdqe,
            qe_ * (1 - se_) / 0.3 - se_ / 3, -qs_ / tdqe,
            qs_ * (1 - ss_) / 0.3 - ss_ / 3,
            dvi, (h_i_inf(vi_) - hi_) / tau_h_i(vi_), (n_i_inf(vi_) - ni_) / tau_n_i(vi_),
            (1 + np.tanh(vi_ / 10)) / 2 * (1 - qi_) / 0.1 - qi_ / tdqi,
            qi_ * (1 - si_) / 0.3 - si_ / 9,
        )

    for k in range(steps):
        old_ve, old_vi = ve.copy(), vi.copy()
        phase = k * dt / pulse_period - pulse_phase if pulse_period else 0.0
        first = slopes(ve, he, ne, qe, se, q_stoch, s_stoch, vi, hi, ni, qi, si, phase)
        mid = tuple(x + 0.5 * dt * dx for x, dx in zip(
            (ve, he, ne, qe, se, q_stoch, s_stoch, vi, hi, ni, qi, si), first
        ))
        phase_mid = (k + 0.5) * dt / pulse_period - pulse_phase if pulse_period else 0.0
        second = slopes(*mid, phase_mid)
        values = tuple(x + dt * dx for x, dx in zip(
            (ve, he, ne, qe, se, q_stoch, s_stoch, vi, hi, ni, qi, si), second
        ))
        ve, he, ne, qe, se, q_stoch, s_stoch, vi, hi, ni, qi, si = values

        crossed_e = np.flatnonzero((old_ve > -20) & (ve <= -20))
        crossed_i = np.flatnonzero((old_vi > -20) & (vi <= -20))
        if crossed_e.size:
            e_indices.extend(crossed_e.tolist())
            e_times.extend((((-20 - ve[crossed_e]) * k * dt + (old_ve[crossed_e] + 20) * (k + 1) * dt)
                            / (old_ve[crossed_e] - ve[crossed_e])).tolist())
        if crossed_i.size:
            i_indices.extend(crossed_i.tolist())
            i_times.extend((((-20 - vi[crossed_i]) * k * dt + (old_vi[crossed_i] + 20) * (k + 1) * dt)
                            / (old_vi[crossed_i] - vi[crossed_i])).tolist())
        q_stoch[rng.random(num_e) < 40 / 1000 * dt] = 1

    return {
        "t_e_spikes": np.asarray(e_times), "i_e_spikes": np.asarray(e_indices, dtype=int),
        "t_i_spikes": np.asarray(i_times), "i_i_spikes": np.asarray(i_indices, dtype=int),
        "f_hat_e": len(e_times) / num_e / t_final * 1000,
        "f_hat_i": len(i_times) / num_i / t_final * 1000,
    }
