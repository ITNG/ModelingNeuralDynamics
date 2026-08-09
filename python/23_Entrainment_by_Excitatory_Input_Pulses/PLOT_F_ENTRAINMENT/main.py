import numpy as np
import matplotlib.pyplot as plt

T = 25.
tau = 20.
epsilon = 0.2


def F(alpha):
    return (alpha + epsilon) * np.exp(-T / tau)


if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot([0, 1 - epsilon], [F(0), F(1 - epsilon)], '-k', linewidth=4)
    plt.plot([1 - epsilon, 1], [0, 0], '-k', linewidth=6)
    plt.plot([1 - epsilon, 1 - epsilon], [0, F(1 - epsilon)], ':k', linewidth=1)
    plt.plot([0, 1], [0, 1], ':k', linewidth=2)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.gca().set_box_aspect(1)
    plt.xticks([0, 1])
    plt.xlabel(r'$\alpha$')
    plt.ylabel('$F$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
