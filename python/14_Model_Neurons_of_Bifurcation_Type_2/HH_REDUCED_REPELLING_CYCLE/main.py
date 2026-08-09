import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

c = 1.
g_k = 36.
g_na = 120.
g_l = 0.3
v_k = -82.
v_na = 45.
v_l = -59.
i_ext = 5.5
dt = 0.01


def alpha_m(v):
    with np.errstate(divide='ignore', invalid='ignore'):
        out = (v + 45) / 10.0 / (1 - exp(-(v + 45) / 10))
    return out if abs(v + 45) > 1e-8 else 1.0


def alpha_n(v):
    return 0.01 * (-60.0 - v) / (exp((-60 - v) / 10) - 1)


def beta_m(v):
    return 4 * exp(-(v + 70) / 18)


def beta_n(v):
    return 0.125 * exp(-(v + 70) / 80)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def f(v):
    '''dv/dt at i_ext=0 of the reduced (m=m_inf(v), h=0.83-n_inf(v)) HH model'''
    return g_na * m_inf(v) ** 3 * (0.83 - n_inf(v)) * (v_na - v) \
        + g_k * n_inf(v) ** 4 * (v_k - v) + g_l * (v_l - v)


def find_fixed_point():
    v_left, v_right = -100., 50.
    while v_right - v_left > 1e-10:
        v_c = (v_left + v_right) / 2
        if (f(v_c) + i_ext) * (f(v_left) + i_ext) > 0:
            v_left = v_c
        else:
            v_right = v_c
    return (v_left + v_right) / 2


def derivative(v, n):
    m = m_inf(v)
    h = 0.83 - n
    dv = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v) + i_ext) / c
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    return dv, dn


def simulate(v0, n0, t_final, reverse=False):
    sign = -1 if reverse else 1
    m_steps = round(t_final / dt)
    v, n = np.zeros(m_steps + 1), np.zeros(m_steps + 1)
    v[0], n[0] = v0, n0
    for k in range(m_steps):
        v_inc, n_inc = derivative(v[k], n[k])
        v_tmp = v[k] + sign * dt / 2 * v_inc
        n_tmp = n[k] + sign * dt / 2 * n_inc
        v_inc, n_inc = derivative(v_tmp, n_tmp)
        v[k + 1] = v[k] + sign * dt * v_inc
        n[k + 1] = n[k] + sign * dt * n_inc
    return v, n


v_c = find_fixed_point()
n_c = n_inf(v_c)

# settled attracting limit cycle: 1000ms forward from v=20, keep last 20ms
v_attr, n_attr = simulate(20., n_c, t_final=1000.)
tail = round(20 / dt)
v_attr, n_attr = v_attr[-tail:], n_attr[-tail:]

# repelling cycle traced backward from just next to the fixed point,
# keep the last 15ms
v_rep, n_rep = simulate(v_c + 0.001, n_c, t_final=1000., reverse=True)
tail = round(15 / dt)
v_rep, n_rep = v_rep[-tail:], n_rep[-tail:]

if __name__ == "__main__":

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    axes[0].plot(v_attr, n_attr, color='k', linewidth=2)
    axes[0].plot(v_rep, n_rep, ':r', linewidth=2)
    axes[0].set_xlim(-100, 50)
    axes[0].set_ylim(0.3, 0.8)
    axes[0].set_box_aspect(1)
    axes[0].set_xlabel('$v$')
    axes[0].set_ylabel('$n$')

    axes[1].plot(v_attr, n_attr, color='k', linewidth=2)
    axes[1].plot(v_c, n_c, '.', markersize=15)
    axes[1].plot(v_rep, n_rep, ':r', linewidth=2)
    axes[1].set_xlim(-75, -55)
    axes[1].set_ylim(0.35, 0.45)
    axes[1].set_box_aspect(1)
    axes[1].set_xlabel('$v$')
    axes[1].set_ylabel('$n$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
