import numpy as np
import matplotlib.pyplot as plt

tau_m = 10.
I = 0.12
g_I = 0.05
tau_I = 5.
dt = 0.01
x = np.arange(1, 10000) / 10000.


def psif(tau_m, I, g_I, tau_I, dt, x):
    '''x can be a 1-d array, in which case so is the output.'''
    v_1 = 0.
    v_2 = x.copy()
    s = 1.
    N = len(x)
    out = np.zeros(N)
    done = np.zeros(N, dtype=bool)

    dt05 = dt / 2

    while not done.all():
        v_1_old = v_1
        v_2_old = v_2.copy()
        v_1_inc = -v_1 / tau_m + I - g_I * s * v_1
        v_2_inc = -v_2 / tau_m + I - g_I * s * v_2
        v_1_tmp = v_1 + dt05 * v_1_inc
        v_2_tmp = v_2 + dt05 * v_2_inc
        s_tmp = s * np.exp(-dt05 / tau_I)
        v_1_inc = -v_1_tmp / tau_m + I - g_I * s_tmp * v_1_tmp
        v_2_inc = -v_2_tmp / tau_m + I - g_I * s_tmp * v_2_tmp
        v_1 = v_1 + dt * v_1_inc
        v_2 = v_2 + dt * v_2_inc

        ind = np.where((v_2 > 1) & ~done)[0]
        out[ind] = (v_1_old * (v_2[ind] - 1) + v_1 * (1 - v_2_old[ind])) / (v_2[ind] - v_2_old[ind])
        done[ind] = True
        s = s * np.exp(-dt / tau_I)

    return out


psi_vec = psif(tau_m, I, g_I, tau_I, dt, x)
max_psi = psi_vec.max()
min_psi = psi_vec.min()

if __name__ == "__main__":
    print(f"max_psi = {max_psi}")
    print(f"min_psi = {min_psi}")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(x, psi_vec, '-k', linewidth=4)
    ax.plot([0, 1], [0, 1], '--k', linewidth=2)
    ax.axis([0, 1, 0, 1])
    ax.set_xlabel('$x$')
    ax.set_ylabel(r'$\psi(x)$')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
