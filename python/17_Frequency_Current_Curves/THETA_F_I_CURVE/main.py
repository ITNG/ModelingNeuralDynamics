import numpy as np
import matplotlib.pyplot as plt

N = 200
I = 1 / 2 + np.arange(N + 1) / N * 0.1

f = 1000 * np.sqrt(2 * I - 1) / np.pi

if __name__ == "__main__":

    plt.figure(figsize=(7, 3.5))
    plt.plot(I, f, '-k', linewidth=2)
    plt.plot([0, 1 / 2], [0, 0], '-k', linewidth=4)
    plt.xlabel('$I$')
    plt.ylabel('$f$')
    plt.xlim(0.4, 0.6)
    plt.ylim(0, 150)
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
