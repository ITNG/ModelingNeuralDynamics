import numpy as np
import matplotlib.pyplot as plt

phi_B_0 = 0.4
dt = 0.01
dt05 = dt / 2


def g(phi):
    p = np.mod(phi, 1)
    return p ** 2 * (1 - p)


def ceiling(x):
    '''like ceil, but rounds an exact integer up to the next integer
    instead of leaving it unchanged (matches matlab's self-defined
    "ceiling", used because phi==an integer means "just spiked", not
    "about to spike").'''
    c = np.ceil(x)
    if c == x:
        c += 1
    return c


def simulate_event_driven(epsilon, phi_B_0=phi_B_0):
    phi_A, phi_B = 0., phi_B_0
    t_vec = [0.]
    psi_vec = [phi_B - phi_A]

    t_final = 6 / epsilon
    t = 0.

    while t < t_final:
        delta_A = ceiling(phi_A) - phi_A
        delta_B = ceiling(phi_B) - phi_B
        if delta_B < delta_A:
            t = t + delta_B
            phi_A = phi_A + delta_B
            phi_A = phi_A + epsilon * g(phi_A)
            phi_B = ceiling(phi_B)
        else:
            t = t + delta_A
            phi_B = phi_B + delta_A
            phi_B = phi_B + epsilon * g(phi_B)
            phi_A = ceiling(phi_A)
        t_vec.append(t)
        psi_vec.append(phi_B - phi_A)

    return np.array(t_vec), np.array(psi_vec)


def simulate_de(epsilon, phi_B_0=phi_B_0):
    t_final = 6 / epsilon
    m_steps = round(t_final / dt)
    psi = np.zeros(m_steps + 1)
    psi[0] = phi_B_0
    for k in range(m_steps):
        psi_inc = epsilon * (g(psi[k]) - g(-psi[k]))
        psi_tmp = psi[k] + dt05 * psi_inc
        psi_inc = epsilon * (g(psi_tmp) - g(-psi_tmp))
        psi[k + 1] = psi[k] + dt * psi_inc
    t = np.arange(m_steps + 1) * dt
    return t, psi


t_vec_1, psi_vec_1 = simulate_event_driven(epsilon=0.5)
t_de_1, psi_de_1 = simulate_de(epsilon=0.5)

t_vec_2, psi_vec_2 = simulate_event_driven(epsilon=0.1)
t_de_2, psi_de_2 = simulate_de(epsilon=0.1)

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    for n in range(len(t_vec_1) - 1):
        axes[0].plot([t_vec_1[n], t_vec_1[n + 1]], [psi_vec_1[n], psi_vec_1[n]], '-k', linewidth=2)
    axes[0].plot(t_de_1, psi_de_1, '-r', linewidth=2)
    axes[0].set_ylabel(r'$\psi$')
    axes[0].set_title(r'$\epsilon=0.5$')
    M = np.ceil(max(psi_vec_1.max(), psi_de_1.max()))
    m = np.floor(min(psi_vec_1.min(), psi_de_1.min()))
    axes[0].axis([0, 6 / 0.5, m, M])

    for n in range(len(t_vec_2) - 1):
        axes[1].plot([t_vec_2[n], t_vec_2[n + 1]], [psi_vec_2[n], psi_vec_2[n]], '-k', linewidth=2)
    axes[1].plot(t_de_2, psi_de_2, '-r', linewidth=2)
    axes[1].set_ylabel(r'$\psi$')
    axes[1].set_xlabel('$t$ [units of $T$]')
    axes[1].set_title(r'$\epsilon=0.1$')
    M = np.ceil(max(psi_vec_2.max(), psi_de_2.max()))
    m = np.floor(min(psi_vec_2.min(), psi_de_2.min()))
    axes[1].axis([0, 6 / 0.1, m, M])

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
