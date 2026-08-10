from math import exp
import numpy as np

c = 1.
g_k, g_na, g_l = 80., 100., 0.1
v_k, v_na, v_l = -100., 50., -67.

epsilon = 1e-5  # relative tolerance for the bisection searches
dt = 0.01
dt05 = dt / 2
t_final = 5000.
m_steps = round(t_final / dt)

v_rev = -75.  # reversal potential of tonic/periodic inhibition


def g(t):
    return np.exp(np.cos(np.pi * t / 25) ** 4) - 1


def m_inf(v):
    alpha_m = 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))
    beta_m = 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)
    return alpha_m / (alpha_m + beta_m)


def alpha_h(v):
    return 0.128 * exp(-(v + 50) / 18)


def beta_h(v):
    return 4. / (1 + exp(-(v + 27) / 5))


def alpha_n(v):
    return 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))


def beta_n(v):
    return 0.5 * exp(-(v + 57) / 40)


def firing_rate(i_ext, g_store):
    '''run 5000 ms of the RTM neuron under drive i_ext and the
    precomputed inhibitory conductance trace g_store, and return the
    resulting firing rate.'''
    v, m, h, n = -70., m_inf(-70.), 0.7, 0.6
    num_spikes = 0
    for k in range(m_steps):
        v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v)
                  + g_store[k] * (v_rev - v) + i_ext) / c
        h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
        n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n

        v_tmp = v + dt05 * v_inc
        m_tmp = m_inf(v_tmp)
        h_tmp = h + dt05 * h_inc
        n_tmp = n + dt05 * n_inc

        v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
                  + g_l * (v_l - v_tmp) + (g_store[k] + g_store[k + 1]) * 0.5 * (v_rev - v_tmp) + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v_new = v + dt * v_inc
        m = m_inf(v_new)
        h = h + dt * h_inc
        n = n + dt * n_inc

        if v_new < -20 and v >= -20:
            num_spikes += 1
        v = v_new

    return num_spikes / t_final * 1000


def bisect_threshold(g_store, target):
    '''bisection search for the drive at which the firing rate first
    reaches (>=) target, i.e. the boundary between f==0 (or f<target)
    and f>=target.'''
    i_ext_left, i_ext_right = 0., 3.
    while (i_ext_right - i_ext_left) / ((i_ext_right + i_ext_left) / 2) > epsilon:
        i_ext = (i_ext_left + i_ext_right) / 2
        f = firing_rate(i_ext, g_store)
        if f < target:
            i_ext_left = i_ext
        else:
            i_ext_right = i_ext
    return (i_ext_left + i_ext_right) / 2


g_bar_vec = 0.15 + np.arange(5) * 0.05
w = np.zeros(5)
t = np.arange(m_steps + 1) * dt

for ijk, g_bar in enumerate(g_bar_vec):
    g_mean = np.mean(g(np.arange(1, 10001) / 10000 * 25))
    g_factor = g_bar / g_mean
    g_store = g(t) * g_factor

    # I_L: onset of firing (f crosses from 0 to positive)
    I_L = bisect_threshold(g_store, target=1e-9)
    # I_R: onset of full-rate (39 Hz) locked firing
    I_R = bisect_threshold(g_store, target=39.)
    w[ijk] = I_R - I_L

ratios = w[:4] / w[1:]

if __name__ == "__main__":
    for g_bar, w_ijk in zip(g_bar_vec, w):
        print(f"g_bar={g_bar:.2f}  w={w_ijk:.6f}")
    print("w =", w)
    print("ratios =", ratios)
