import numpy as np
import matplotlib.pyplot as plt

tau_m = 10.
I = 0.12
g_I = 0.015
tau_I = 5.
dt = 0.01
x = np.arange(1, 10000) / 10000.


def psif(tau_m, I, g_I, tau_I, dt, x):
    '''x can be a 1-d array, in which case so is the output. Unlike
    PLOT_PHI/PLOT_PSI's psif, the synaptic drive s(t) is a fixed alpha
    function of absolute time (not reset/renormalized to 1 at t=0 with
    exponential decay tied to spike count), matching the matlab source.'''
    v_1 = 0.
    v_2 = x.copy()
    N = len(x)
    out = np.zeros(N)
    done = np.zeros(N, dtype=bool)

    delta = 0.

    def alpha(t):
        return np.exp(-(t - delta) / tau_I) * (t >= delta)

    dt05 = dt / 2
    t = 0.
    s = alpha(0.)

    while not done.all():
        v_1_old = v_1
        v_2_old = v_2.copy()
        v_1_inc = -v_1 / tau_m + I - g_I * s * v_1
        v_2_inc = -v_2 / tau_m + I - g_I * s * v_2
        v_1_tmp = v_1 + dt05 * v_1_inc
        v_2_tmp = v_2 + dt05 * v_2_inc
        s_tmp = alpha(t + dt05)
        v_1_inc = -v_1_tmp / tau_m + I - g_I * s_tmp * v_1_tmp
        v_2_inc = -v_2_tmp / tau_m + I - g_I * s_tmp * v_2_tmp
        v_1 = v_1 + dt * v_1_inc
        v_2 = v_2 + dt * v_2_inc

        ind = np.where((v_2 > 1) & ~done)[0]
        out[ind] = (v_1_old * (v_2[ind] - 1) + v_1 * (1 - v_2_old[ind])) / (v_2[ind] - v_2_old[ind])
        done[ind] = True
        t = t + dt
        s = alpha(t)

    return out


psi_vec = psif(tau_m, I, g_I, tau_I, dt, x)
phi_vec = psif(tau_m, I, g_I, tau_I, dt, psi_vec)

if __name__ == "__main__":
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    ax = axes[0]
    ax.plot(x, psi_vec, '.k', markersize=2)
    ax.plot([0, 1], [0, 1], '--k', linewidth=1)
    ax.axis([0, 1, 0, 1])
    ax.set_xlabel('$x$')
    ax.set_ylabel(r'$\psi$')
    ax.set_aspect('equal')

    ax = axes[1]
    ax.plot(x, phi_vec, '.k', markersize=2)
    ax.plot([0, 1], [0, 1], '--k', linewidth=1)
    ax.axis([0, 1, 0, 1])
    ax.set_xlabel('$x$')
    ax.set_ylabel(r'$\phi$')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
