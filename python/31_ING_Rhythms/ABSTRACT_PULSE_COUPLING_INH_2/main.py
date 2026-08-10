import numpy as np
import matplotlib.pyplot as plt

a = 4.
epsilon = 0.4


def g(phi):
    return -epsilon * phi * np.tanh((1 - phi) * a) / np.tanh(a)


varphi = np.arange(1001) / 1000

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(varphi, g(varphi), '-k', linewidth=5)
    plt.axis([0, 1, -1, 0])
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\varphi$')
    plt.ylabel(r'$g(\varphi)$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
