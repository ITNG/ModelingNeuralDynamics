import numpy as np
import pylab as pl
from numpy import exp


def a_inf(v):
    return 1.0 / (1.0 + exp(-(v + 14.0) / 16.6))


def b_inf(v):
    return 1.0 / (1.0 + exp((v + 71.0) / 7.3))


def tau_a(v):
    return [5.0] * len(v)


def tau_b(v):
    return 1.0 / (0.000009 / exp((v - 26.0) / 28.5) +
                  0.014/(0.2+exp(-(v+70.0)/11.0)))


if __name__ == "__main__":

    fig, ax = pl.subplots(nrows=2, ncols=2,
                          figsize=(6, 3), sharex=True)
    v = np.arange(-100, 50, 0.01)

    ax[0, 0].plot(v, a_inf(v), color='k', lw=2)
    ax[0, 0].set_ylabel(r"$a_{\infty}(v)$")

    ax[0, 1].plot(v, tau_a(v), color='k', lw=2)
    ax[0, 1].set_ylabel(r"$\tau_{a}$ [ms]")

    ax[1, 0].plot(v, b_inf(v), color='k', lw=2)
    ax[1, 0].set_ylabel(r"$b_{\infty}(v)$")

    ax[1, 1].plot(v, tau_b(v), color='k', lw=2)
    ax[1, 1].set_ylabel(r"$\tau_{b}$ [ms]")
    # ax[1, 1].set_ylim(0, 10)

    for i in range(2):
        for j in range(2):
            ax[i, j].tick_params(labelsize=10)

    for i in range(2):
        ax[i, 0].set_ylim(0, 1)
        ax[i, 0].set_yticks([0, 0.5, 1])

    ax[0, 0].set_xlim(np.min(v), np.max(v))

    for i in range(2):
        ax[1, i].set_xlabel("v [mV]")
        ax[1, i].set_xticks([-100, -50, 0, 50])

    pl.tight_layout()
    pl.savefig("fig_34_7.png")
    pl.close()
    # pl.show()
