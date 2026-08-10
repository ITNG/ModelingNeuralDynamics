import numpy as np
import matplotlib.pyplot as plt

t = np.arange(-250, 251) / 100.


def phi_of(alpha):
    phi = np.exp(alpha * np.sin(np.pi * t) ** 2) - 1
    return phi / phi[1:].mean()


if __name__ == "__main__":
    fig, axes = plt.subplots(3, 1, figsize=(6, 8))

    # alpha=1e-5 approximates alpha->0 without dividing by zero
    for ax, alpha, label in zip(axes, [1e-5, 5, 10], [r'$\alpha \rightarrow 0$', r'$\alpha=5$', r'$\alpha=10$']):
        ax.plot(t, phi_of(alpha), '-k', linewidth=2)
        ax.set_title(label)
        ax.axis([-2.5, 2.5, 0, 7])

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
