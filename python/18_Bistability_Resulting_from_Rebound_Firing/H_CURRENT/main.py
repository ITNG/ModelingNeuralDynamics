import numpy as np
import matplotlib.pyplot as plt


def r_inf(v):
    return 1. / (1 + np.exp((v + 84) / 10.2))


def tau_r(v):
    return 1. / (np.exp(-14.59 - 0.086 * v) + np.exp(-1.87 + 0.0701 * v))


v = np.arange(-100, 51)

if __name__ == "__main__":

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(v, r_inf(v), '-k', linewidth=2)
    ax[0].set_xlabel('$v$ [mV]')
    ax[0].set_ylabel(r'$r_\infty(v)$')
    ax[0].set_xlim(-100, 50)
    ax[0].set_ylim(0, 1)

    ax[1].plot(v, tau_r(v), '-k', linewidth=2)
    ax[1].set_xlabel('$v$ [mV]')
    ax[1].set_ylabel(r'$\tau_r$ [ms]')
    ax[1].set_xlim(-100, 50)
    ax[1].set_ylim(0, 1000)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
