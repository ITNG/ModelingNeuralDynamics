from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

a = 1.25
tau_n = 15.625
i_ext = -0.5
t_final = 400.
dt = 0.01


def derivative(x0, t):
    '''
    define the FitzHugh-Nagumo model
    '''
    v, n = x0
    dv = v - v ** 3 / 3 - n + i_ext
    dn = (a * v - n) / tau_n
    return [dv, dn]


x0 = [-1., -2.]

if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v, n = sol[:, 0], sol[:, 1]

    v_line = np.arange(-100, 101) / 100. * 3
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

    ax[0].plot(v_line, a * v_line, color='k', linewidth=2)
    ax[0].plot(v_line, v_line - v_line ** 3 / 3 + i_ext, color='r', linewidth=2)
    ax[0].plot(v, n, color='b', linewidth=2)
    ax[0].set_xlim(-3, 3)
    ax[0].set_ylim(-3, 3)
    ax[0].set_xlabel('$v$')
    ax[0].set_ylabel('$n$')

    ax[1].plot(t, v, color='k', linewidth=2)
    ax[1].set_xlim(0, t_final)
    ax[1].set_ylim(-3, 3)
    ax[1].set_xlabel('$t$')
    ax[1].set_ylabel('$v$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
