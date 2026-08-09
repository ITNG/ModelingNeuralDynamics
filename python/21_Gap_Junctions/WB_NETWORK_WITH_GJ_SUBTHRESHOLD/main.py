import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 9.
g_na = 35.
g_l = 0.1
v_k = -90.
v_na = 55.
v_l = -65.

num = 2
gap_strength = 0.01

i_ext = np.array([1., 0.])

t_final = 100.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


def alpha_h(v):
    return 0.35 * exp(-(v + 58) / 20)


def alpha_m(v):
    return 0.1 * (v + 35) / (1 - exp(-(v + 35) / 10))


def alpha_n(v):
    return 0.05 * (v + 34) / (1 - exp(-0.1 * (v + 34)))


def beta_h(v):
    return 5. / (exp(-0.1 * (v + 28)) + 1)


def beta_m(v):
    return 4 * exp(-(v + 60) / 18)


def beta_n(v):
    return 0.625 * exp(-(v + 44) / 80)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def simulate(gap_gate):
    '''forward-Euler (Heun) integration of two WB neurons coupled by a gap
    junction. gap_gate(v1) returns the 2x2 coupling matrix as a function
    of neuron 1's current voltage -- constant for the always-on gap
    junction, zeroed out above threshold for the subthreshold-only one.'''
    v = np.zeros((num, m_steps + 1))
    m = np.zeros((num, m_steps + 1))
    h = np.zeros((num, m_steps + 1))
    n = np.zeros((num, m_steps + 1))
    v[:, 0] = [-63., -63.]
    m[:, 0] = m_inf(v[:, 0])
    h[:, 0] = h_inf(v[:, 0])
    n[:, 0] = n_inf(v[:, 0])

    for k in range(m_steps):
        g_gap = gap_gate(v[0, k])
        d_gap = g_gap.sum(axis=0)

        v_inc = (g_k * n[:, k] ** 4 * (v_k - v[:, k]) + g_na * m[:, k] ** 3 * h[:, k] * (v_na - v[:, k])
                  + g_l * (v_l - v[:, k]) + g_gap @ v[:, k] - d_gap * v[:, k] + i_ext) / c
        n_inc = alpha_n(v[:, k]) * (1 - n[:, k]) - beta_n(v[:, k]) * n[:, k]
        h_inc = alpha_h(v[:, k]) * (1 - h[:, k]) - beta_h(v[:, k]) * h[:, k]

        v_tmp = v[:, k] + dt05 * v_inc
        h_tmp = h[:, k] + dt05 * h_inc
        n_tmp = n[:, k] + dt05 * n_inc
        m_tmp = m_inf(v_tmp)

        v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
                  + g_l * (v_l - v_tmp) + g_gap @ v_tmp - d_gap * v_tmp + i_ext) / c
        h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
        n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

        v[:, k + 1] = v[:, k] + dt * v_inc
        h[:, k + 1] = h[:, k] + dt * h_inc
        n[:, k + 1] = n[:, k] + dt * n_inc
        m[:, k + 1] = m_inf(v[:, k + 1])

    return v


gap_always_on = np.array([[0., gap_strength], [gap_strength, 0.]])
v_always = simulate(lambda v1: gap_always_on)


def gap_subthreshold(v1, thr=-50.):
    return gap_always_on if v1 <= thr else np.zeros((2, 2))


v_sub = simulate(gap_subthreshold)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt
    half = m_steps // 2

    fig, ax = plt.subplots(2, figsize=(7, 6))

    ax[0].plot(t, v_always[0], '-k', linewidth=3)
    ax[0].plot(t, v_sub[0], '-r', linewidth=1)
    ax[0].set_ylabel('$v_1$ [mV]')
    ax[0].set_xlim(0, t_final)
    ax[0].set_ylim(-100, 50)

    ax[1].plot(t, v_always[1], '-k', linewidth=3)
    ax[1].plot(t, v_sub[1], '-r', linewidth=1)
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$v_2$ [mV]')
    ax[1].set_xlim(0, t_final)
    tail = v_always[1, half:]
    ax[1].set_ylim(tail.min() - 0.75, tail.max() + 0.75)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
