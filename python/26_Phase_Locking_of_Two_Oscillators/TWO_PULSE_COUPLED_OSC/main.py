import numpy as np
import matplotlib.pyplot as plt


def g(phi):
    return phi ** 2 * (1 - phi)


def f(phi):
    return phi + g(phi)


phi_A = 0.
phi_B = 0.5
num_spikes_A = 1
num_spikes_B = 0
t_spikes_A = [0.]
t_spikes_B = []

N = 12
t = 0.

for k in range(N):
    t = t + (1 - phi_B)
    num_spikes_B += 1
    t_spikes_B.append(t)
    phi_A = f(1 - phi_B)
    phi_B = 0.
    t = t + (1 - phi_A)
    num_spikes_A += 1
    t_spikes_A.append(t)
    phi_B = f(1 - phi_A)
    phi_A = 0.

t_spikes_A = np.array(t_spikes_A)
t_spikes_B = np.array(t_spikes_B)

if __name__ == "__main__":

    plt.figure(figsize=(8, 4))
    plt.plot(t_spikes_A, np.ones(num_spikes_A), '.r', markersize=15)
    plt.plot(t_spikes_B, 2 * np.ones(num_spikes_B), '.b', markersize=15)
    plt.axis([0, N - 2, 0, 3])
    plt.xticks(range(1, N))
    plt.yticks([])
    plt.xlabel('$t$ [units of $T$]')
    plt.title('spike times of A (red) and B (blue)')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
