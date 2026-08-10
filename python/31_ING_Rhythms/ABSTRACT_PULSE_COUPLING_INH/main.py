import numpy as np
import matplotlib.pyplot as plt


def H(s):
    return (1 + np.tanh(s)) / 2


def g(phi):
    return -0.5 * phi * (H((-phi + 0.8) / 0.05) - H(-4.)) * H((phi - 0.1) / 0.1)


def f(phi):
    return phi + g(phi)


def bigF(phi):
    return f(1 - phi)


def bigG(phi):
    return bigF(bigF(phi))


phi = np.arange(1001) / 1000

phi_left = phi[:1000]
phi_right = phi[1:1001]
G_left = bigG(phi_left)
G_right = bigG(phi_right)
ind = np.where((G_left - phi_left) * (G_right - phi_right) <= 0)[0]

if __name__ == "__main__":

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))

    axes[0].plot(phi, g(phi), '-k', linewidth=2)
    axes[0].axis([0, 1, -1, 0])
    axes[0].set_box_aspect(1)
    axes[0].set_xlabel(r'$\varphi$')
    axes[0].set_ylabel('$g$')

    ax = axes[1]
    ax.plot(phi, bigG(phi), '-k', linewidth=2)
    ax.plot([0, 1], [0, 1], '--k', linewidth=1)
    ax.plot(phi_left[ind[0]], G_left[ind[0]], '.g', markersize=15)
    ax.plot(phi_right[ind[4]], G_right[ind[4]], '.g', markersize=15)
    ax.plot(phi_left[ind[1]], G_left[ind[1]], 'or', markersize=8, linewidth=2, fillstyle='none')
    ax.plot(phi_right[ind[3]], G_right[ind[3]], 'or', markersize=8, linewidth=2, fillstyle='none')
    ax.plot((phi_left[ind[2]] + phi_right[ind[2]]) / 2, (G_left[ind[2]] + G_right[ind[2]]) / 2,
            '.b', markersize=15)
    ax.axis([0, 1, 0, 1])
    ax.set_box_aspect(1)
    ax.set_xlabel(r'$\varphi$')
    ax.set_ylabel('$G$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
