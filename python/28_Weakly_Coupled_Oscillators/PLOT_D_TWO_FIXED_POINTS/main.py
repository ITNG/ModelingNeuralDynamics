import numpy as np
import matplotlib.pyplot as plt


def g_0(phi):
    phi_tilde = np.mod(phi, 1)
    return phi_tilde ** 2 * (1 - phi_tilde)


def D(psi):
    return g_0(psi) - g_0(1 - psi)


psi = np.arange(301) / 300
c = 0.08

psi_0 = 0.5
for i in range(1, 101):
    psi_i = 0.5 + 0.5 * i / 100
    if D(psi_i) > c:
        psi_0 = psi_i

psi_left, psi_right = 0.5, psi_0
while psi_right - psi_left > 1e-12:
    psi_c = (psi_left + psi_right) / 2
    if D(psi_c) > c:
        psi_right = psi_c
    else:
        psi_left = psi_c
psi_unstable = (psi_left + psi_right) / 2

psi_left, psi_right = psi_0, 1.
while psi_right - psi_left > 1e-12:
    psi_c = (psi_left + psi_right) / 2
    if D(psi_c) > c:
        psi_left = psi_c
    else:
        psi_right = psi_c
psi_stable = (psi_left + psi_right) / 2

if __name__ == "__main__":

    plt.figure(figsize=(6, 6))
    plt.plot(psi, D(psi), '-k', linewidth=5)
    plt.axis([0, 1, -0.1, 0.1])
    plt.gca().set_box_aspect(1)
    plt.xlabel(r'$\psi$')
    plt.ylabel(r'$D(\psi)$')

    plt.plot([0, 1], [c, c], '--r', linewidth=3)
    plt.text(-0.075, 0.08, '$c$', color='red', fontsize=20)

    plt.plot([psi_unstable, psi_unstable], [-0.1, D(psi_unstable)], '--r', linewidth=2)
    plt.plot(psi_unstable, -0.1, 'or', markersize=12, markerfacecolor='w')

    plt.plot([psi_stable, psi_stable], [-0.1, D(psi_stable)], '--r', linewidth=2)
    plt.plot(psi_stable, -0.1, 'or', markersize=12, markerfacecolor='r')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
