import numpy as np
import matplotlib.pyplot as plt


def g(phi):
    return phi ** 2 * (1 - phi)


def f(phi):
    return phi + g(phi)


def bigF(phi):
    return f(1 - phi)


def bigG(phi):
    return bigF(bigF(phi))


phi = np.arange(101) / 100

f_left = f(phi[:100])
f_right = f(phi[1:101])
if np.min(f_right - f_left) <= 0:
    print('f is not strictly increasing')

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 2, figsize=(8, 8))

    axes[0, 0].plot(phi, g(phi), '-k', linewidth=2)
    axes[0, 0].axis([0, 1, 0, 1])
    axes[0, 0].set_box_aspect(1)
    axes[0, 0].set_xlabel(r'$\varphi$')
    axes[0, 0].set_ylabel('$g$')

    axes[0, 1].plot(phi, f(phi), '-k', linewidth=2)
    axes[0, 1].axis([0, 1, 0, 1])
    axes[0, 1].set_box_aspect(1)
    axes[0, 1].set_xlabel(r'$\varphi$')
    axes[0, 1].set_ylabel('$f$')

    axes[1, 0].plot(phi, bigF(phi), '-k', linewidth=2)
    axes[1, 0].axis([0, 1, 0, 1])
    axes[1, 0].set_box_aspect(1)
    axes[1, 0].set_xlabel(r'$\varphi$')
    axes[1, 0].set_ylabel('$F$')

    axes[1, 1].plot(phi, bigG(phi), '-k', linewidth=2)
    axes[1, 1].plot([0, 1], [0, 1], '--k', linewidth=1)
    axes[1, 1].axis([0, 1, 0, 1])
    axes[1, 1].set_box_aspect(1)
    axes[1, 1].set_xlabel(r'$\varphi$')
    axes[1, 1].set_ylabel('$G$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
