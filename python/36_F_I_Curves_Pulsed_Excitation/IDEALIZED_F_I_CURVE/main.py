import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(6, 3))

    I_below = np.arange(-100, 1) / 100
    ax.plot(I_below, np.zeros_like(I_below), '-k', linewidth=4)

    I_above = np.arange(0, 101) / 100
    ax.plot(I_above, np.sqrt(I_above), '-k', linewidth=2)

    ax.set_ylabel('$f$')
    ax.text(-0.02, -0.1, '$I_c$')
    ax.axis([-1, 1, 0, 1])
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
