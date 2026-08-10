import numpy as np
import matplotlib.pyplot as plt

alpha = 1.
Period = 25.
g_bar = 0.1
N = 2000
m = np.mean(np.exp(alpha * np.cos(np.pi * np.arange(N) / N) ** 2) - 1)
tau = 10.
I_vec = np.arange(1, 101) / 100 * 0.2 + 0.1

t_final = 2000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)


def g(t):
    return g_bar * (np.exp(alpha * np.cos(np.pi * t / Period) ** 2) - 1) / m


f_vec = np.zeros(len(I_vec))
for ij, I in enumerate(I_vec):
    v = 0.
    num_spikes = 0
    for k in range(1, m_steps + 1):
        v_inc = -v / tau + I - g((k - 1) * dt) * v
        v_tmp = v + dt05 * v_inc
        v_inc = -v_tmp / tau + I - g((k - 0.5) * dt) * v_tmp
        v = v + dt * v_inc
        if v > 1:
            v = 0.
            num_spikes += 1
    f_vec[ij] = num_spikes / t_final * 1000

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(7, 5))

    L = len(I_vec)
    for k in range(L - 1):
        if f_vec[k + 1] - f_vec[k] <= 1:
            ax.plot([I_vec[k], I_vec[k + 1]], [f_vec[k], f_vec[k + 1]], '-b', linewidth=2)
        else:
            mid = (I_vec[k] + I_vec[k + 1]) / 2
            ax.plot([mid, mid], [f_vec[k], f_vec[k + 1]], '--b', linewidth=1)

    # closed-form F-I curve of the LIF neuron under a constant (mean)
    # inhibitory conductance g_bar
    f = np.arange(1, 201)
    A = f / 1000 / (1 + tau * g_bar)
    A = 1 / A
    A = A / tau
    A = np.exp(A)
    I_closed = A * (1 + tau * g_bar) / tau / (A - 1)
    ax.plot(I_closed, f, '-r', linewidth=2)
    ax.axis([I_vec.min(), I_vec.max(), 0, 200])
    ax.set_xlabel('$I$')
    ax.set_ylabel('$f$ [Hz]')

    I_onset_tonic = 1 / tau + g_bar
    ax.plot([I_vec.min(), I_onset_tonic], [0, 0], '-r', linewidth=4)
    ax.plot(I_onset_tonic, 0, '.r', markersize=20)

    for k in range(L - 1):
        if f_vec[k] == 0 and f_vec[k + 1] > 0:
            I_onset = (I_vec[k] + I_vec[k + 1]) / 2
            ax.plot(I_onset, 0, '.b', markersize=20)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
