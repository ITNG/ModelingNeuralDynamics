import numpy as np
import matplotlib.pyplot as plt

w_IE = 1.
w_EI = 1.
w_II = 0.
tau_E = 5.
tau_I = 10.
I_E = 20.
I_I = 0.

t_final = 1000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


def f(x):
    return 100 * x ** 2 / (900 + x ** 2) * (x > 0)


def g(x):
    return 100 * x ** 2 / (400 + x ** 2) * (x > 0)


def simulate(w_EE):
    E = np.zeros(m_steps + 1)
    I = np.zeros(m_steps + 1)
    E[0], I[0] = 10., 10.
    for k in range(m_steps):
        E_inc = (f(w_EE * E[k] - w_IE * I[k] + I_E) - E[k]) / tau_E
        I_inc = (g(w_EI * E[k] - w_II * I[k] + I_I) - I[k]) / tau_I
        E_tmp = E[k] + dt05 * E_inc
        I_tmp = I[k] + dt05 * I_inc
        E_inc = (f(w_EE * E_tmp - w_IE * I_tmp + I_E) - E_tmp) / tau_E
        I_inc = (g(w_EI * E_tmp - w_II * I_tmp + I_I) - I_tmp) / tau_I
        E[k + 1] = E[k] + dt * E_inc
        I[k + 1] = I[k] + dt * I_inc
    return E, I


ind = slice(round(4 * m_steps / 5) - 1, m_steps + 1)  # matlab's ind is 1-indexed

w_EE_vec = [1.5, 1.25, 1.0]
panels = []
for w_EE in w_EE_vec:
    E, I = simulate(w_EE)
    panels.append((w_EE, E[ind], I[ind]))

if __name__ == "__main__":

    fig, ax = plt.subplots(1, 3, figsize=(12, 4.5))

    for a, (w_EE, E, I) in zip(ax, panels):
        a.plot(E, I, '-k', linewidth=2)
        a.set_xlim(0, 100)
        a.set_ylim(0, 100)
        a.set_box_aspect(1)
        a.set_xlabel('$E$')
        a.set_title(f'$w_{{EE}}={w_EE}$')

    ax[0].set_ylabel('$I$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
