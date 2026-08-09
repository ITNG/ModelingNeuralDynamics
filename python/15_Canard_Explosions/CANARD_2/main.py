import numpy as np
import matplotlib.pyplot as plt

a = 5.
tau_n = 60.
I = -4.256889

t_final = 10000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

v = np.zeros(m_steps + 1)
n = np.zeros(m_steps + 1)
v[0], n[0] = 0., 0.

for k in range(m_steps):
    v_inc = v[k] - v[k] ** 3 / 3 - n[k] + I
    n_inc = (a * v[k] - n[k]) / tau_n
    v_tmp = v[k] + dt05 * v_inc
    n_tmp = n[k] + dt05 * n_inc
    v_inc = v_tmp - v_tmp ** 3 / 3 - n_tmp + I
    n_inc = (a * v_tmp - n_tmp) / tau_n
    v[k + 1] = v[k] + dt * v_inc
    n[k + 1] = n[k] + dt * n_inc

tail_start = 4 * m_steps // 5 - 1  # matlab's v(4*m_steps/5 : m_steps+1) is 1-indexed
vv = v[tail_start:]
nn = n[tail_start:]

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(vv, nn, '-b', linewidth=3)

    A, B, C, D = -2, 0.5, -5, -4
    v_plot = A + np.arange(101) / 100 * (B - A)
    plt.plot(v_plot, v_plot - v_plot ** 3 / 3 + I, '-g', linewidth=2)
    plt.plot(v_plot, a * v_plot, '-r', linewidth=2)

    plt.xlim(A, B)
    plt.ylim(C, D)
    plt.gca().set_box_aspect(1)
    plt.xlabel('$v$')
    plt.ylabel('$n$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
