import numpy as np
import matplotlib.pyplot as plt

tau = 10.
num = 2
g_gap = 0.01
beta = 3.

G = np.array([[0., g_gap], [g_gap, 0.]])
c = G.sum(axis=0)

i_ext = np.array([0.125, 0.09])

t_final = 100.
dt = 0.002
m_steps = round(t_final / dt)


def simulate(epsilon):
    '''forward-Euler integration of two LIF neurons coupled by a gap
    junction (diffusive term G*v-c*v), plus an extra epsilon kick to the
    other cell's voltage whenever a cell spikes (epsilon=0 isolates the
    pure diffusive gap-junction coupling from that spike-triggered kick)'''
    v = np.zeros((num, m_steps + 1))
    v[:, 0] = [0.4, 0.9]
    spike_times = [[], []]

    for k in range(m_steps):
        v_inc = -v[:, k] / tau + i_ext + G @ v[:, k] - c * v[:, k]
        v[:, k + 1] = v[:, k] + dt * v_inc
        w = v[:, k + 1]
        ind = int(np.argmax(w))
        if w[ind] > 1:
            v[ind, k + 1] = 0.
            other = 1 - ind
            v[other, k + 1] += epsilon
            if v[other, k + 1] > 1:
                v[other, k + 1] = 0.
            spike_times[ind].append((k + 1) * dt)

    return v, spike_times


v_coupled, spikes_coupled = simulate(beta * g_gap)
v_uncoupled, spikes_uncoupled = simulate(0.)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    fig, ax = plt.subplots(2, figsize=(7, 6))

    for tt in spikes_coupled[0]:
        ax[0].plot([tt, tt], [0, 6], '-k', linewidth=3)
    for tt in spikes_uncoupled[0]:
        ax[0].plot([tt, tt], [0, 6], '-r', linewidth=1)
    ax[0].plot(t, v_coupled[0], '-k', linewidth=3)
    ax[0].plot(t, v_uncoupled[0], '-r', linewidth=1)
    ax[0].set_ylabel('$v_1$ [mV]')
    ax[0].set_xlim(0, t_final)
    ax[0].set_ylim(0, 6)

    for tt in spikes_coupled[1]:
        ax[1].plot([tt, tt], [0, 6], '-k', linewidth=3)
    for tt in spikes_uncoupled[1]:
        ax[1].plot([tt, tt], [0, 6], '-r', linewidth=1)
    ax[1].plot(t, v_coupled[1], '-k', linewidth=3)
    ax[1].plot(t, v_uncoupled[1], '-r', linewidth=1)
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$v_2$ [mV]')
    ax[1].set_xlim(0, t_final)
    ax[1].set_ylim(0.85, 0.95)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
