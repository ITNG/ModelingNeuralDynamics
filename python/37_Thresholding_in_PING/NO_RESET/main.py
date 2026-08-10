import numpy as np
import matplotlib.pyplot as plt

tau_m = 10.
I = 0.12

T = 25.
epsilon = 10.
g_bar = 1. / 7.
tau_m_hat = 1. / (1. / tau_m + g_bar * T / epsilon)

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(8, 4))

    v0 = 1.0
    for ijk in range(5):
        t = np.arange(101) / 100 * epsilon
        v = v0 * np.exp(-t / tau_m_hat) + tau_m_hat * I * (1 - np.exp(-t / tau_m_hat))
        ax.plot(t + ijk * T, v, '-k', linewidth=2)
        v0 = v[-1]

        t = np.arange(101) / 100 * (T - epsilon)
        v = v0 * np.exp(-t / tau_m) + tau_m * I * (1 - np.exp(-t / tau_m))
        ax.plot(t + ijk * T + epsilon, v, '-k', linewidth=2)
        v0 = v[-1]

    ax.set_xlabel('$t$')
    ax.set_ylabel('$v$')
    ax.axis([0, 5 * T, 0, 1])

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
