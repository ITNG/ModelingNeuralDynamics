import numpy as np
import matplotlib.pyplot as plt

N = 200
phi = np.arange(N + 1) / N
tau_m = 2.
I = 0.13

g_hat = 1 / np.pi / np.sqrt(tau_m * I - 1 / 4) / (1 + np.tan(np.pi * (phi - 1 / 2)) ** 2)

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(phi, g_hat, '-k', linewidth=6)
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\varphi$')
    plt.ylabel(r'$\hat{g}$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
