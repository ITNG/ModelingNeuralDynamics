import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 80.
g_na = 100.
g_l = 0.1
v_k = -100.
v_na = 50.
v_l = -67.

i_ext = 1.2
N = 10
t_final = 30.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


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


def rtm_init(i_ext, phi_vec):
    '''find the RTM limit cycle at i_ext (plain-float Heun integration
    until the 5th spike), then interpolate (v, h, n) at phases phi_vec
    (fraction of the last full period, measured from the 4th spike).'''
    t_final_init = 5000.
    dt_init = 0.01
    dt05_init = dt_init / 2

    v = [-70.]
    m = [m_inf(v[0])]
    h = [h_inf(v[0])]
    n = [n_inf(v[0])]
    t_spikes = []

    k = 0
    t = 0.
    while len(t_spikes) < 5 and t < t_final_init:
        vk, mk, hk, nk = v[k], m[k], h[k], n[k]
        v_inc = (g_k * nk ** 4 * (v_k - vk) + g_na * mk ** 3 * hk * (v_na - vk) + g_l * (v_l - vk) + i_ext) / c
        h_inc = alpha_h(vk) * (1 - hk) - beta_h(vk) * hk
        n_inc = alpha_n(vk) * (1 - nk) - beta_n(vk) * nk

        v_tmp = vk + dt05_init * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = hk + dt05_init * h_inc
        n_tmp = nk + dt05_init * n_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v.append(vk + dt_init * v_inc)
        h.append(hk + dt_init * h_inc)
        n.append(nk + dt_init * n_inc)
        m.append(m_inf(v[-1]))

        if vk >= -20 and v[-1] < -20:
            t_spike = (k * dt_init * (-20 - v[-1]) + (k + 1) * dt_init * (20 + vk)) / (vk - v[-1])
            t_spikes.append(t_spike)
        t = (k + 1) * dt_init
        k += 1

    num = len(phi_vec)
    T = t_spikes[4] - t_spikes[3]
    out = np.zeros((num, 3))
    for i, phi0 in enumerate(phi_vec):
        t0 = phi0 * T + t_spikes[3]
        kk = int(t0 / dt_init)
        frac_hi = (t0 - kk * dt_init) / dt_init
        frac_lo = ((kk + 1) * dt_init - t0) / dt_init
        out[i, 0] = v[kk + 1] * frac_hi + v[kk] * frac_lo
        out[i, 1] = h[kk + 1] * frac_hi + h[kk] * frac_lo
        out[i, 2] = n[kk + 1] * frac_hi + n[kk] * frac_lo
    return out


phi_vec = np.arange(N - 1, -1, -1) / N + 1 / (2 * N)


def simulate(g_syn, v_syn, tau_syn):
    initial_vector = rtm_init(i_ext, phi_vec)
    v = initial_vector[:, 0].copy()
    m = m_inf(v)
    h = initial_vector[:, 1].copy()
    n = initial_vector[:, 2].copy()

    v_trace = np.zeros((N, m_steps + 1))
    v_trace[:, 0] = v

    for k in range(m_steps):
        t = k * dt
        v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v)
                  + i_ext + g_syn * exp(-t / tau_syn) * (v_syn - v)) / c
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h

        v_tmp = v + dt05 * v_inc
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc
        m_tmp = m_inf(v_tmp)

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext + g_syn * exp(-(t + dt05) / tau_syn) * (v_syn - v_tmp)) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v = v + dt * v_inc
        h = h + dt * h_inc
        n = n + dt * n_inc
        m = m_inf(v)

        v_trace[:, k + 1] = v

    return v_trace


v_trace_none = simulate(g_syn=0., v_syn=-75., tau_syn=1.)  # g_syn=0 -> no inhibition
v_trace_weak = simulate(g_syn=0.30, v_syn=-75., tau_syn=9.)
v_trace_strong = simulate(g_syn=2.25, v_syn=-75., tau_syn=1.)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt
    fig, axes = plt.subplots(3, 1, figsize=(8, 9))

    for ax, v_trace in zip(axes, (v_trace_none, v_trace_weak, v_trace_strong)):
        for k in range(N):
            ax.plot(t, v_trace[k], '-k', linewidth=1)
        ax.axis([0, t_final, -100, 50])
        ax.set_ylabel('$v$ [mV]')

    axes[-1].set_xlabel('$t$ [ms]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
