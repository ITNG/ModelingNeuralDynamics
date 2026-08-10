import numpy as np
import matplotlib.pyplot as plt

dt = 0.01
dt05 = dt / 2


def P0(tau_m, J, g, tau_I):
    I = J + 1 / tau_m
    v = 0.
    k = 1
    while v < 1:
        v_old = v
        v_inc = -v / tau_m + I - g * np.exp(-(k - 1) * dt / tau_I) * v
        v_tmp = v + dt05 * v_inc
        v_inc = -v_tmp / tau_m + I - g * np.exp(-(k - 0.5) * dt / tau_I) * v_tmp
        v = v + dt * v_inc
        k += 1
    period = ((k - 2) * dt * (v - 1) + (k - 1) * dt * (1 - v_old)) / (v - v_old)
    return period


def P1(tau_m, J, g, tau_I):
    I = J + 1 / tau_m
    if g <= J:
        return 0.
    v = 1.
    k = 1
    while v <= 1:
        v_old = v
        v_inc = -v / tau_m + I - g * np.exp(-(k - 1) * dt / tau_I) * v
        v_tmp = v + dt05 * v_inc
        v_inc = -v_tmp / tau_m + I - g * np.exp(-(k - 0.5) * dt / tau_I) * v_tmp
        v = v + dt * v_inc
        k += 1
    period = ((k - 2) * dt * (v - 1) + (k - 1) * dt * (1 - v_old)) / (v - v_old)
    return period


def P(tau_m, J, g, tau_I):
    return (P0(tau_m, J, g, tau_I) + P1(tau_m, J, g, tau_I)) / 2


def S(tau_m, J, g, tau_I):
    return (P0(tau_m, J, g, tau_I) - P1(tau_m, J, g, tau_I)) / P(tau_m, J, g, tau_I)


tau_m_0, J_0, g_0, tau_I_0 = 10., 0.02, 0.15, 9.

J_vec = J_0 * 0.8 + np.arange(101) / 100 * (J_0 * 1.2 - J_0 * 0.8)
tau_I_vec = tau_I_0 * 0.8 + np.arange(101) / 100 * tau_I_0 * 0.4
g_vec = g_0 * 0.8 + np.arange(101) / 100 * g_0 * 0.4

P_vec_tau_I = np.array([P(tau_m_0, J_0, g_0, tI) for tI in tau_I_vec])
S_vec_tau_I = np.array([S(tau_m_0, J_0, g_0, tI) for tI in tau_I_vec])

P_vec_g = np.array([P(tau_m_0, J_0, gg, tau_I_0) for gg in g_vec])
S_vec_g = np.array([S(tau_m_0, J_0, gg, tau_I_0) for gg in g_vec])

P_vec_J = np.array([P(tau_m_0, JJ, g_0, tau_I_0) for JJ in J_vec])
S_vec_J = np.array([S(tau_m_0, JJ, g_0, tau_I_0) for JJ in J_vec])

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))

    axes[0, 0].plot(tau_I_vec, P_vec_tau_I, '-k', linewidth=2)
    axes[0, 0].plot(tau_I_0, P(tau_m_0, J_0, g_0, tau_I_0), '.r', markersize=15)
    axes[0, 0].axis([tau_I_0 * 0.8, tau_I_0 * 1.2, 20, 40])
    axes[0, 0].set_ylabel('$P$')
    axes[0, 0].set_title(r'varying $\tau_I$')

    axes[1, 0].plot(tau_I_vec, S_vec_tau_I, '-k', linewidth=2)
    axes[1, 0].plot(tau_I_0, S(tau_m_0, J_0, g_0, tau_I_0), '.r', markersize=15)
    axes[1, 0].axis([tau_I_0 * 0.8, tau_I_0 * 1.2, 0, 0.1])
    axes[1, 0].set_xlabel(r'$\tau_I$')
    axes[1, 0].set_ylabel('$S$')

    axes[0, 1].plot(g_vec, P_vec_g, '-k', linewidth=2)
    axes[0, 1].plot(g_0, P(tau_m_0, J_0, g_0, tau_I_0), '.r', markersize=15)
    axes[0, 1].set_xticks(np.arange(0.13, 0.171, 0.02))
    axes[0, 1].axis([0.8 * g_0, 1.2 * g_0, 20, 40])
    axes[0, 1].set_title(r'varying $\overline{g}_{\rm syn}$')

    axes[1, 1].plot(g_vec, S_vec_g, '-k', linewidth=2)
    axes[1, 1].plot(g_0, S(tau_m_0, J_0, g_0, tau_I_0), '.r', markersize=15)
    axes[1, 1].set_xticks(np.arange(0.13, 0.171, 0.02))
    axes[1, 1].set_xlabel(r'$\overline{g}_{\rm syn}$')
    axes[1, 1].axis([0.8 * g_0, 1.2 * g_0, 0, 0.1])

    axes[0, 2].plot(J_vec, P_vec_J, '-k', linewidth=2)
    axes[0, 2].plot(J_0, P(tau_m_0, J_0, g_0, tau_I_0), '.r', markersize=15)
    axes[0, 2].set_xticks([0.017, 0.02, 0.023])
    axes[0, 2].axis([0.8 * J_0, 1.2 * J_0, 20, 40])
    axes[0, 2].set_title('varying $J$')

    axes[1, 2].plot(J_vec, S_vec_J, '-k', linewidth=2)
    axes[1, 2].plot(J_0, S(tau_m_0, J_0, g_0, tau_I_0), '.r', markersize=15)
    axes[1, 2].set_xticks([0.017, 0.02, 0.023])
    axes[1, 2].set_xlabel('$J$')
    axes[1, 2].axis([0.8 * J_0, 1.2 * J_0, 0, 0.1])

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
