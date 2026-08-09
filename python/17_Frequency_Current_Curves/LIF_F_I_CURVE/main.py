import numpy as np
import matplotlib.pyplot as plt

N = 1000
ind = np.arange(1, N + 1)
I = 0.1 + np.exp(-N / ind) * np.exp(1) * 0.1

tau_m = 10.

T = tau_m * np.log(tau_m * I / (tau_m * I - 1))
f = 1000. / T

if __name__ == "__main__":

    plt.figure(figsize=(7, 3.5))
    plt.plot(I, f, '-k', linewidth=2)
    plt.plot([0, 0.1], [0, 0], '-k', linewidth=4)
    plt.xlabel('$I$')
    plt.ylabel('$f$')
    plt.xlim(0, 0.2)
    plt.ylim(0, 150)
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
