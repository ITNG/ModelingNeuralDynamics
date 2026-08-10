import numpy as np
from numpy import exp, tanh
import matplotlib.pyplot as plt

c = 1.
g_k, g_na, g_l = 9., 35., 0.1
v_k, v_na, v_l = -90., 55., -65.

U = 0.5
C = 1.25

tau_rec = 500.
tau_d_q = 5.
tau_r = 3.
tau_d = 9.

i_ext = 0.5

t_final = 200.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


# Note: unlike the WB neuron used elsewhere in this project (which
# divides tau_h/tau_n by a phi=5 speedup factor), this particular
# matlab source (alpha_h.m/beta_h.m/alpha_n.m/beta_n.m) bakes phi=5
# directly into the rate constants themselves (e.g. alpha_h=0.35*exp(...)
# = 5 * 0.07*exp(...)). alpha_m/beta_m are unaffected, as usual.


def alpha_m(v):
    return 0.1 * (v + 35) / (1 - exp(-(v + 35) / 10))


def beta_m(v):
    return 4. * exp(-(v + 60) / 18)


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def alpha_h(v):
    return 0.35 * exp(-(v + 58) / 20)


def beta_h(v):
    return 5. / (exp(-0.1 * (v + 28)) + 1)


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def alpha_n(v):
    return 0.05 * (v + 34) / (1 - exp(-0.1 * (v + 34)))


def beta_n(v):
    return 0.625 * exp(-(v + 44) / 80)


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


v = np.zeros(m_steps + 1)
m = np.zeros(m_steps + 1)
h = np.zeros(m_steps + 1)
n = np.zeros(m_steps + 1)
p = np.zeros(m_steps + 1)
q = np.zeros(m_steps + 1)
s = np.zeros(m_steps + 1)

v[0] = -70.
m[0] = m_inf(v[0])
h[0] = h_inf(v[0])
n[0] = n_inf(v[0])
p[0] = 1.
q[0] = 0.
s[0] = 0.

num_spikes = 0

for k in range(m_steps):
    v_inc = (g_k * n[k] ** 4 * (v_k - v[k]) + g_na * m[k] ** 3 * h[k] * (v_na - v[k])
              + g_l * (v_l - v[k]) + i_ext) / c
    n_inc = alpha_n(v[k]) * (1 - n[k]) - beta_n(v[k]) * n[k]
    h_inc = alpha_h(v[k]) * (1 - h[k]) - beta_h(v[k]) * h[k]
    p_inc = -C * (1 + tanh(v[k] / 10)) * p[k] * np.log(1 / (1 - U)) + (1 - p[k] - q[k]) / tau_rec
    q_inc = C * (1 + tanh(v[k] / 10)) * p[k] * np.log(1 / (1 - U)) - q[k] / tau_d_q
    s_inc = q[k] * (1 - s[k]) / tau_r - s[k] / tau_d

    v_tmp = v[k] + dt05 * v_inc
    h_tmp = h[k] + dt05 * h_inc
    n_tmp = n[k] + dt05 * n_inc
    m_tmp = m_inf(v_tmp)
    p_tmp = p[k] + dt05 * p_inc
    q_tmp = q[k] + dt05 * q_inc
    s_tmp = s[k] + dt05 * s_inc

    v_inc = (g_k * n_tmp ** 4 * (v_k - v_tmp) + g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp)
              + g_l * (v_l - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp
    p_inc = -C * (1 + tanh(v_tmp / 10)) * p_tmp * np.log(1 / (1 - U)) + (1 - p_tmp - q_tmp) / tau_rec
    q_inc = C * (1 + tanh(v_tmp / 10)) * p_tmp * np.log(1 / (1 - U)) - q_tmp / tau_d_q
    s_inc = q_tmp * (1 - s_tmp) / tau_r - s_tmp / tau_d

    v[k + 1] = v[k] + dt * v_inc
    h[k + 1] = h[k] + dt * h_inc
    n[k + 1] = n[k] + dt * n_inc
    m[k + 1] = m_inf(v[k + 1])
    p[k + 1] = p[k] + dt * p_inc
    q[k + 1] = q[k] + dt * q_inc
    s[k + 1] = s[k] + dt * s_inc

    if v[k + 1] < -20 and v[k] >= -20:
        num_spikes += 1

gamma = C * (1 + tanh(v / 10))
integral_of_delta_function_per_period = gamma.sum() * dt / num_spikes

if __name__ == "__main__":
    print(f"integral_of_delta_function_per_period = {integral_of_delta_function_per_period}")

    t = np.arange(m_steps + 1) * dt

    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    axes[0, 0].plot(t, v, '-k', linewidth=2)
    axes[0, 0].set_ylabel('$v$ [mV]')
    axes[0, 0].axis([0, t_final, -100, 50])

    axes[0, 1].plot(t, p, '-k', linewidth=2)
    axes[0, 1].set_ylabel('$p$')
    axes[0, 1].axis([0, t_final, 0, 1])

    axes[1, 0].plot(t, q, '-k', linewidth=2)
    axes[1, 0].set_xlabel('$t$ [ms]')
    axes[1, 0].set_ylabel('$q$')
    axes[1, 0].axis([0, t_final, 0, q.max() * 1.1])

    axes[1, 1].plot(t, s, '-k', linewidth=2)
    axes[1, 1].set_xlabel('$t$ [ms]')
    axes[1, 1].set_ylabel('$s$')
    axes[1, 1].axis([0, t_final, 0, s.max() * 1.1])

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
