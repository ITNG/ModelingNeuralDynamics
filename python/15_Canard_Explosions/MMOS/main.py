import numpy as np
import matplotlib.pyplot as plt

a = 5.
tau_n = 60.
I_ext = -4.2
tau_adapt = 150.
delta = 0.2

t_final = 1000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

v = np.zeros(m_steps + 1)
n = np.zeros(m_steps + 1)
v[0], n[0] = -1., -4.75
I_adapt = -delta

for k in range(m_steps):
    v_inc = v[k] - v[k] ** 3 / 3 - n[k] + I_ext + I_adapt
    n_inc = (a * v[k] - n[k]) / tau_n
    v_tmp = v[k] + dt05 * v_inc
    n_tmp = n[k] + dt05 * n_inc
    v_inc = v_tmp - v_tmp ** 3 / 3 - n_tmp + I_ext + I_adapt * np.exp(-dt05 / tau_adapt)
    n_inc = (a * v_tmp - n_tmp) / tau_n
    v[k + 1] = v[k] + dt * v_inc
    n[k + 1] = n[k] + dt * n_inc
    I_adapt = I_adapt * np.exp(-dt / tau_adapt)
    if v[k + 1] < 0 and v[k] >= 0:
        I_adapt = I_adapt - delta

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    plt.figure(figsize=(7, 3.5))
    plt.plot(t, v, '-k', linewidth=2)
    plt.xlim(0, t_final)
    plt.ylim(-3, 3)
    plt.xlabel('$t$')
    plt.ylabel('$v$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
