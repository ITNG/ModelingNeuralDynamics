import numpy as np
import matplotlib.pyplot as plt

epsilon = 1 / 3


def g(phi):
    return epsilon * phi * (1 - phi)


def f(phi):
    return phi + g(phi)


def simulate(delta, t_final=20.):
    '''event-driven simulation of two pulse-coupled oscillators A, B with
    a transmission delay delta (in units of the intrinsic period): each
    spike's effect on the other oscillator arrives delta time units
    later, rather than instantaneously.'''
    phi_A, phi_B = 0., 0.9
    t_A_to_B = delta  # time until the next A-input reaches B
    t_B_to_A = np.inf  # time until the next B-input reaches A
    t_present = 0.

    num_spikes_A, num_spikes_B = 1, 0
    t_spikes_A = [0.]
    t_spikes_B = []

    while t_present < t_final:
        T_vec = [1 - phi_A, 1 - phi_B, t_A_to_B, t_B_to_A]
        T_0 = min(T_vec)
        done = False

        if T_0 == 1 - phi_A:
            phi_B = phi_B + 1 - phi_A
            t_B_to_A = t_B_to_A - (1 - phi_A)
            t_A_to_B = delta
            t_present = t_present + 1 - phi_A
            num_spikes_A += 1
            t_spikes_A.append(t_present)
            phi_A = 0.
            done = True

        if T_0 == 1 - phi_B and not done:
            phi_A = phi_A + 1 - phi_B
            t_A_to_B = t_A_to_B - (1 - phi_B)
            t_B_to_A = delta
            t_present = t_present + 1 - phi_B
            num_spikes_B += 1
            t_spikes_B.append(t_present)
            phi_B = 0.
            done = True

        if T_0 == t_A_to_B and not done:
            phi_B = f(phi_B + t_A_to_B)
            phi_A = phi_A + t_A_to_B
            t_B_to_A = t_B_to_A - t_A_to_B
            t_present = t_present + t_A_to_B
            t_A_to_B = np.inf
            done = True

        if T_0 == t_B_to_A and not done:
            phi_A = f(phi_A + t_B_to_A)
            phi_B = phi_B + t_B_to_A
            t_A_to_B = t_A_to_B - t_B_to_A
            t_present = t_present + t_B_to_A
            t_B_to_A = np.inf
            done = True

    return np.array(t_spikes_A), np.array(t_spikes_B)


t_final = 20.
t_spikes_A_1, t_spikes_B_1 = simulate(delta=0.1, t_final=t_final)
t_spikes_A_2, t_spikes_B_2 = simulate(delta=0.7, t_final=t_final)

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 1, figsize=(8, 6))

    axes[0].plot(t_spikes_A_1, np.ones(len(t_spikes_A_1)), '.r', markersize=15)
    axes[0].plot(t_spikes_B_1, 2 * np.ones(len(t_spikes_B_1)), '.g', markersize=15)
    axes[0].axis([t_final - 10, t_final, 0, 3])
    axes[0].set_yticks([])
    axes[0].set_title(r'$\delta=0.1$:')

    axes[1].plot(t_spikes_A_2, np.ones(len(t_spikes_A_2)), '.r', markersize=15)
    axes[1].plot(t_spikes_B_2, 2 * np.ones(len(t_spikes_B_2)), '.g', markersize=15)
    axes[1].axis([t_final - 10, t_final, 0, 3])
    axes[1].set_yticks([])
    axes[1].set_xlabel('$t$ [units of $T$]')
    axes[1].set_title(r'$\delta=0.7$:')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
