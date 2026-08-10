import numpy as np
import matplotlib.pyplot as plt


def f_tilde(s):
    return 0.5 * (s + 1 / np.sqrt(2)) * (1 / np.sqrt(2) - s)


s = np.zeros(101)
phi = np.zeros(101)
f = np.zeros(101)
for k in range(101):
    s[k] = -1 / np.sqrt(2) + k / 100 * np.sqrt(2)
    u = f_tilde(s[k]) / np.sqrt(2)
    phi[k] = 0.5 + s[k] / np.sqrt(2) - u
    f[k] = 0.5 + s[k] / np.sqrt(2) + u

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(phi, f, '-k', linewidth=6)
    plt.axis([0, 1, 0, 1])
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\varphi$')
    plt.ylabel('$f$')
    plt.plot([0, 1], [0, 1], '--k', linewidth=2)

    epsilon = 0.025
    plt.plot([0.5 - epsilon, 0.5 + epsilon], [0.5 + epsilon, 0.5 - epsilon], '-k', linewidth=1)
    plt.plot([0.7 - epsilon, 0.7 + epsilon], [0.7 + epsilon, 0.7 - epsilon], '-k', linewidth=1)
    plt.text(0.55, 0.43, '$0$', fontsize=20, rotation=45)
    plt.text(0.75, 0.63, '$s$', fontsize=20, rotation=45)
    delta = 0.15
    plt.plot([0.7, 0.7 - delta], [0.7, 0.7 + delta], '-r', linewidth=4)
    plt.plot([0.7, 0.7 - delta], [0.7, 0.7], '-b', linewidth=2)
    plt.plot([0.7 - delta, 0.7 - delta], [0.7, 0.7 + delta], '-b', linewidth=2)
    plt.text(0.49, 0.75, '$u$', color='b', fontsize=20)
    plt.text(0.59, 0.67, '$u$', color='b', fontsize=20)
    plt.text(0.67, 0.8, r'$\tilde{f}$', color='r', fontsize=20)
    plt.plot(0.7 - delta, 0.7 + delta, '.k', markersize=20)
    plt.text(0.7 - delta - 0.43, 0.7 + delta + 0.05, r'$(\varphi,f(\varphi))$', fontsize=20)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
