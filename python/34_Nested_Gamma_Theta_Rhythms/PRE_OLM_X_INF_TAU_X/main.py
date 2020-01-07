import numpy as np
import pylab as pl

import numpy as np
import pylab as pl
from numpy import exp


def alpha_h(v):
    return 0.07 * exp(-(v + 63.0) / 20.0)


def alpha_m(v):
    q = (v+38.0)/10.0
    return q / (1.0 - exp(-q))


def alpha_n(v):
    return 0.018 * (v - 25.0) / (1.0 - exp(-(v - 25.0) / 25.0))


def beta_h(v):
    return 1.0 / (exp(-(v + 33.0) / 10.0) + 1.0)


def beta_m(v):
    return 4.0 * exp(-(v + 65.0) / 18.0)


def beta_n(v):
    return 0.0036 * (35.0 - v) / (1.0 - exp(-(35.0 - v) / 12.0))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def tau_h(v):
    return 1.0 / (alpha_h(v) + beta_h(v))


def tau_m(v):
    return 1.0 / (alpha_m(v) + beta_m(v))


def tau_n(v):
    return 1.0 / (alpha_n(v) + beta_n(v))


if __name__ == "__main__":


    fig, ax = pl.subplots(nrows=3, ncols=2,
    figsize=(6, 5.5), sharex=True)
    v = np.arange(-100, 50, 0.01)

    ax[0, 0].plot(v, m_inf(v), color='k', lw=2)
    ax[0, 0].set_ylabel(r"$m_{\infty}$")

    # ax[0, 1].plot(v, tau_m(v), color='k', lw=2)
    # ax[0, 1].set_ylabel(r"$\tau_{m}$ [ms]")
    ax[0, 1].axis("off")

    ax[1, 0].plot(v, h_inf(v), color='k', lw=2)
    ax[1, 0].set_ylabel(r"$h_{\infty}$")

    ax[1, 1].plot(v, tau_h(v), color='k', lw=2)
    ax[1, 1].set_ylabel(r"$\tau_{h}$ [ms]")
    ax[1, 1].set_ylim(0, 10)

    ax[2, 0].plot(v, n_inf(v), color='k', lw=2)
    ax[2, 0].set_ylabel(r"$n_{\infty}$")

    ax[2, 1].plot(v, tau_n(v), color='k', lw=2)
    ax[2, 1].set_ylabel(r"$\tau_{n}$ [ms]")
    ax[2, 1].set_ylim(0, 4)



    for i in range(3):
        for j in range(2):
            ax[i, j].tick_params(labelsize=10)

    for i in range(3):
        ax[i, 0].set_ylim(0, 1)
        ax[i, 0].set_yticks([0, 0.5, 1])

    ax[0, 0].set_xlim(np.min(v), np.max(v))
    
    for i in range(2):
        ax[2, i].set_xlabel("v [mV]")
        ax[2, i].set_xticks([-100, -50, 0, 50])


    pl.tight_layout()
    pl.savefig("fig_34_4.png")
    pl.close()
    # pl.show()
