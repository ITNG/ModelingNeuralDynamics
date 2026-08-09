import numpy as np
import matplotlib.pyplot as plt

phi = np.arange(101) / 100
tau_m = 2.
I = 0.13
dv = 0.1

g = np.arctan(np.tan(np.pi * phi - np.pi / 2) + dv / np.sqrt(tau_m * I - 1 / 4)) / np.pi + 1 / 2 - phi

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(phi, g, '-k', linewidth=6)
    plt.axis([0, 1, 0, 1])
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\varphi$')
    plt.ylabel('$g$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
