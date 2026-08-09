import numpy as np
import matplotlib.pyplot as plt

v = np.arange(-100, 51)
B = 1. / (1 + np.exp(-0.062 * v) / 3.57)

if __name__ == "__main__":

    plt.figure(figsize=(7, 3.5))
    plt.plot(v, B)
    plt.xlabel(r'$v_{\rm post}$')
    plt.ylabel('$B$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
