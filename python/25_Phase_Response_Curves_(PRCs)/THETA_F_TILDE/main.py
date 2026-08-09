import numpy as np
import matplotlib.pyplot as plt

phi = np.arange(101) / 100
tau_m = 2.
I = 0.13
dv = 0.1

f = np.arctan(np.tan(np.pi * phi - np.pi / 2) + dv / np.sqrt(tau_m * I - 1 / 4)) / np.pi + 1 / 2

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(phi, f, '-k', linewidth=6)
    plt.plot([0, 1], [0, 1], '--k', linewidth=2)

    epsilon = 0.025
    plt.plot([0.5 - epsilon, 0.5 + epsilon], [0.5 + epsilon, 0.5 - epsilon], '-k', linewidth=1)
    plt.plot([0.7 - epsilon, 0.7 + epsilon], [0.7 + epsilon, 0.7 - epsilon], '-k', linewidth=1)
    plt.text(0.55, 0.43, '$0$', fontsize=20, rotation=45)
    plt.text(0.75, 0.63, '$s$', fontsize=20, rotation=45)
    delta = 0.09
    plt.plot([0.7, 0.7 - delta], [0.7, 0.7 + delta], '-r', linewidth=4)
    plt.text(0.5, 0.62, r'$\tilde{f}(s)$', color='red', rotation=45, fontsize=20)

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\varphi$')
    plt.ylabel('$f$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
