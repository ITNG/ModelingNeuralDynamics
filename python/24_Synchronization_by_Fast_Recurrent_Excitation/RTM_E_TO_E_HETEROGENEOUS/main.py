import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

# matlab's rng('default'); rng(63806) seeds MATLAB's own Mersenne-Twister
# stream, which numpy cannot reproduce bit-for-bit -- this port uses its
# own numpy seed instead, so the specific realization of random i_ext/
# g_syn (and hence the exact spike times) differs from matlab's, even
# though the model and algorithm are identical. Verified by shape/
# statistics, not by an exact numeric match (see the test file).
rng = np.random.default_rng(63806)

c = 1.
g_k = 80.
g_na = 100.
g_l = 0.1
v_k = -100.
v_na = 50.
v_l = -67.

N = 30
i_ext = 0.25 + rng.random(N) * 0.1
g_syn = 0.00625 + rng.random((N, N)) * 0.0025
np.fill_diagonal(g_syn, 0.)  # no autapses

tau_r = 0.5
tau_peak = 0.5
tau_d = 2.


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


def rtm_init(i_ext_i, phi):
    '''find the RTM limit cycle at i_ext_i (plain-float Heun integration
    until the 5th spike), then interpolate (v, h, n) at phase phi
    (fraction of the last full period, measured from the 4th spike).'''
    t_final = 5000.
    dt = 0.001
    dt05 = dt / 2

    v = [-70.]
    m = [m_inf(v[0])]
    h = [h_inf(v[0])]
    n = [n_inf(v[0])]
    t_spikes = []

    k = 0
    t = 0.
    while len(t_spikes) < 5 and t < t_final:
        vk, mk, hk, nk = v[k], m[k], h[k], n[k]
        v_inc = (g_k * nk ** 4 * (v_k - vk) + g_na * mk ** 3 * hk * (v_na - vk)
                  + g_l * (v_l - vk) + i_ext_i) / c
        h_inc = alpha_h(vk) * (1 - hk) - beta_h(vk) * hk
        n_inc = alpha_n(vk) * (1 - nk) - beta_n(vk) * nk

        v_tmp = vk + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = hk + dt05 * h_inc
        n_tmp = nk + dt05 * n_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext_i) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v.append(vk + dt * v_inc)
        h.append(hk + dt * h_inc)
        n.append(nk + dt * n_inc)
        m.append(m_inf(v[-1]))

        if vk >= -20 and v[-1] < -20:
            t_spike = ((k) * dt * (-20 - v[-1]) + (k + 1) * dt * (20 + vk)) / (vk - v[-1])
            t_spikes.append(t_spike)
        t = (k + 1) * dt
        k += 1

    if len(t_spikes) < 5:
        return np.array([v[k], h[k], n[k]])

    T = t_spikes[4] - t_spikes[3]
    t0 = phi * T + t_spikes[3]
    kk = int(t0 / dt)
    frac_hi = (t0 - kk * dt) / dt
    frac_lo = ((kk + 1) * dt - t0) / dt
    return np.array([
        v[kk + 1] * frac_hi + v[kk] * frac_lo,
        h[kk + 1] * frac_hi + h[kk] * frac_lo,
        n[kk + 1] * frac_hi + n[kk] * frac_lo,
    ])


t_final = 200.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

v = np.zeros(N)
h = np.zeros(N)
n = np.zeros(N)
for i in range(N):
    v[i], h[i], n[i] = rtm_init(i_ext[i], rng.random())
m = m_inf(v)
q = np.zeros(N)
s = np.zeros(N)

t_spikes = []
i_spikes = []

g_syn_T = g_syn.T

for k in range(1, m_steps + 1):
    v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v)
              + g_l * (v_l - v) - (g_syn_T @ s) * v + i_ext) / c
    h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
    n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
    q_inc = 5 * (1 + np.tanh(v / 10)) * (1 - q) - q / tau_dq
    s_inc = q * (1 - s) / tau_r - s / tau_d

    v_tmp = v + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h + dt05 * h_inc
    n_tmp = n + dt05 * n_inc
    q_tmp = q + dt05 * q_inc
    s_tmp = s + dt05 * s_inc

    v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
              + g_l * (v_l - v_tmp) - (g_syn_T @ s_tmp) * v_tmp + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
    q_inc = 5 * (1 + np.tanh(v_tmp / 10)) * (1 - q_tmp) - q_tmp / tau_dq
    s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d

    v_old = v.copy()
    v = v + dt * v_inc
    m = m_inf(v)
    h = h + dt * h_inc
    n = n + dt * n_inc
    q = q + dt * q_inc
    s = s + dt * s_inc

    ind = np.where((v_old >= -20) & (v < -20))[0]
    if len(ind) > 0:
        t_new = ((k - 1) * dt * (-v[ind] - 20) + k * dt * (20 + v_old[ind])) / (v_old[ind] - v[ind])
        t_spikes.extend(t_new)
        i_spikes.extend(ind + 1)  # 1-indexed neuron number, matching matlab

t_spikes = np.array(t_spikes)
i_spikes = np.array(i_spikes)

ind1 = np.where(i_spikes == 1)[0]
frequency = 1000 / (t_spikes[ind1[-1]] - t_spikes[ind1[-2]]) if len(ind1) >= 2 else np.nan

if __name__ == "__main__":

    print(f"frequency = {frequency}")

    plt.figure(figsize=(7, 4))
    plt.plot(t_spikes, i_spikes, '.r', markersize=10)
    plt.xlim(t_final - 200, t_final)
    plt.ylim(0, N + 1)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('neuron #')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
