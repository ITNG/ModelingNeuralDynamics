import numpy as np
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

T = 25.
tau = 50.
epsilon = 0.6


def F(alpha):
    return (alpha + epsilon) * np.exp(-T / tau)


alpha = [0.]
alpha.append(F(alpha[0]))
alpha.append(F(alpha[1]))

if __name__ == "__main__":

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.plot([0, 1 - epsilon], [F(0), F(1 - epsilon)], '-k', linewidth=4)
    ax.plot([1 - epsilon, 1], [0, 0], '-k', linewidth=8)
    ax.plot([0, 1], [0, 1], ':k', linewidth=2)

    a1, a2, a3 = alpha
    ax.plot([0, 0], [0, a2], '-r', linewidth=2)
    ax.plot([0, a2], [a2, a2], '-r', linewidth=1)
    ax.plot([a2, a2], [a2, a3], '-r', linewidth=1)
    ax.plot([a2, a3], [a3, a3], '-r', linewidth=1)
    ax.plot([a3, a3], [a3, 0], '-r', linewidth=1)
    ax.plot([0, a3], [0, 0], '-r', linewidth=2)
    ax.plot([a2, a2], [0, a2], ':r', linewidth=1)
    ax.plot(a1, 0, '.r', markersize=20)
    ax.plot(a2, 0, '.r', markersize=20)
    ax.plot(a3, 0, '.r', markersize=20)

    xlim, ylim = (0, 1), (0, 1)
    arrows = [
        (0., a2 / 2, (0., 1.)),
        (a2 / 2, a2, (1., 0.)),
        (a2, a2 + (a3 - a2) / 2, (0., 1.)),
        ((a2 + a3) / 2, a3, (1., 0.)),
        (a3, a3 / 2, (0., -1.)),
        (a3 * 0.4, 0., (-1., 0.)),
    ]
    for x, y, v in arrows:
        draw_arrow(ax, xlim, ylim, x, y, np.array(v), epsilon=0.03, width=2, color='r')

    ax.text(-0.02, -0.08, r'$\alpha_1$', fontsize=18, color='r')
    ax.text(a2 - 0.02, -0.08, r'$\alpha_2$', fontsize=18, color='r')
    ax.text(a3 - 0.02, -0.08, r'$\alpha_3$', fontsize=18, color='r')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1)
    ax.set_xticks([1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel('$F$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
