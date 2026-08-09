import numpy as np
import matplotlib.pyplot as plt


def tau_r_original(v):
    return 1. / (np.exp(-14.59 - 0.086 * v) + np.exp(-1.87 + 0.0701 * v))


def tau_r_modified(v):
    tau_r = tau_r_original(v)
    factor = (0.05 * v + 5) * (v < -80) + (v >= -80)
    return factor ** 2 * tau_r


v = np.arange(-100, 51)

if __name__ == "__main__":

    plt.figure(figsize=(7, 5))
    plt.plot(v, tau_r_modified(v), '--b', linewidth=2)
    plt.plot(v, tau_r_original(v), '-k', linewidth=2)
    plt.xlabel('$v$ [mV]')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
