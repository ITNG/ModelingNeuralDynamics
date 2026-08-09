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

i_ext = 0.30


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
    (fraction of the last full period, measured from the 4th spike).
    Returns (num, 3) array of initial conditions and the period T.'''
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
        v_inc = (g_k * nk ** 4 * (v_k - vk) + g_na * mk ** 3 * hk * (v_na - vk) + g_l * (v_l - vk) + i_ext) / c
        h_inc = alpha_h(vk) * (1 - hk) - beta_h(vk) * hk
        n_inc = alpha_n(vk) * (1 - nk) - beta_n(vk) * nk

        v_tmp = vk + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = hk + dt05 * h_inc
        n_tmp = nk + dt05 * n_inc

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + i_ext) / c
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

    num = len(phi_vec)
    if len(t_spikes) < 5:
        v_last, h_last, n_last = v[k], h[k], n[k]
        out = np.array([[v_last, h_last, n_last]] * num)
        return out, np.inf

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


N = 200
delta_v = 4.
dt = 0.001
dt05 = dt / 2

phi_vec = np.arange(1, N + 1) / N - 1 / (2 * N)

initial_vector, T = rtm_init(i_ext, phi_vec)

v = initial_vector[:, 0].copy() + delta_v  # instantaneous voltage kick (short input pulse)
m = m_inf(v)
h = initial_vector[:, 1].copy()
n = initial_vector[:, 2].copy()

t_star = np.full(N, np.nan)
num_spikes = np.zeros(N, dtype=int)

k = 0
while np.min(num_spikes) < 1:
    k += 1

    v_inc = (g_k * n ** 4 * (v_k - v) + g_na * m ** 3 * h * (v_na - v) + g_l * (v_l - v) + i_ext) / c
    h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
    n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n

    v_tmp = v + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h + dt05 * h_inc
    n_tmp = n + dt05 * n_inc

    v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
              + g_l * (v_l - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

    v_old = v.copy()
    v = v + dt * v_inc
    m = m_inf(v)
    h = h + dt * h_inc
    n = n + dt * n_inc

    ind = np.where((v_old >= -20) & (v < -20))[0]
    for i in ind:
        if num_spikes[i] == 0:
            t_star[i] = ((k - 1) * dt * (-20 - v[i]) + k * dt * (v_old[i] + 20)) / (v_old[i] - v[i])
        num_spikes[i] += 1

g_vec = -t_star / T - phi_vec + 1

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(phi_vec, g_vec, '-k', linewidth=6)
    plt.plot([0, 1], [1, 0], '--k', linewidth=2)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\varphi$')
    plt.ylabel('$g$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
