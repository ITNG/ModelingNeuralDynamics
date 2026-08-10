import numpy as np
import matplotlib.pyplot as plt

epsilon = 2.


def g(phi):
    return epsilon * phi * (1 - phi) ** 3


def f(phi):
    return phi + g(phi)


def bigF(phi):
    return f(1 - phi)


def bigG(phi):
    return bigF(bigF(phi))


phi = np.arange(101) / 100

if __name__ == "__main__":

    fig, axes = plt.subplots(1, 2, figsize=(8, 4.5))

    axes[0].plot(phi, g(phi), '-k', linewidth=2)
    axes[0].axis([0, 1, 0, 1])
    axes[0].set_box_aspect(1)
    axes[0].set_xlabel(r'$\varphi$')
    axes[0].set_ylabel('$g$')

    axes[1].plot(phi, bigG(phi), '-k', linewidth=2)
    axes[1].plot([0, 1], [0, 1], '--k', linewidth=1)
    axes[1].axis([0, 1, 0, 1])
    axes[1].set_box_aspect(1)
    axes[1].set_xlabel(r'$\varphi$')
    axes[1].set_ylabel('$G$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
