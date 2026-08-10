import numpy as np
import matplotlib.pyplot as plt

# MATLAB's rng('default'); rng(63806) cannot be bit-reproduced by NumPy,
# so we use our own seed and verify statistically/visually instead of
# expecting an exact match to MATLAB's spike times.
rng = np.random.default_rng(63806)

N = 3


def g(phi, scale):
    return phi * (1 - phi) / 3 / scale


def simulate(delta, phi0, t_final=200., scale=1.):
    '''event-driven simulation of N all-to-all pulse-coupled oscillators
    with a uniform transmission delay delta (in units of the intrinsic
    period). Each spike sends N-1 signals (one to every other
    oscillator), each arriving delta time units later and nudging the
    receiving oscillator's phase by g(phi)/scale.'''
    t_0 = 0.
    phi = phi0.copy()
    num_spikes = 0
    t_spikes = []
    i_spikes = []

    M = 0
    t = np.zeros(N * N)
    k_to = np.zeros(N * N, dtype=int)
    k_from = np.zeros(N * N, dtype=int)

    while t_0 < t_final:
        next_spike = np.min(1 - phi)
        next_signal = np.inf if M == 0 else np.min(t[:M])
        next_event = min(next_spike, next_signal)

        if next_event == next_signal:
            j0 = np.max(np.where(t[:M] == next_signal)[0])
            phi = phi + next_signal
            t_0 = t_0 + next_signal
            t = t - next_signal
            target = k_to[j0]
            phi[target] = phi[target] + g(phi[target], scale)

            t[j0:M - 1] = t[j0 + 1:M]
            k_to[j0:M - 1] = k_to[j0 + 1:M]
            k_from[j0:M - 1] = k_from[j0 + 1:M]
            M -= 1

        else:
            i0 = np.max(np.where(1 - phi == next_spike)[0])
            phi = phi + next_spike
            t_0 = t_0 + next_spike
            t[:M] = t[:M] - next_spike
            phi[i0] = 0.

            num_spikes += 1
            t_spikes.append(t_0)
            i_spikes.append(i0)

            targets = [j for j in range(N) if j != i0]
            for l, tgt in enumerate(targets):
                k_from[M + l] = i0
                k_to[M + l] = tgt
                t[M + l] = delta
            M += N - 1

    return np.array(t_spikes), np.array(i_spikes)


t_final = 200.
phi0_1 = 0.9 + rng.random(N) * 0.1
t_spikes_1, i_spikes_1 = simulate(delta=0.45, phi0=phi0_1, t_final=t_final, scale=N - 1)

phi0_2 = 0.9 + rng.random(N) * 0.1
t_spikes_2, i_spikes_2 = simulate(delta=0.55, phi0=phi0_2, t_final=t_final, scale=1.)

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 1, figsize=(8, 6))

    axes[0].plot(t_spikes_1, i_spikes_1 + 1, '.k', markersize=15)
    axes[0].axis([t_final - 10, t_final, 0, N + 1])
    axes[0].set_yticks([1, 2, 3])
    axes[0].set_title(r'$\delta=0.45$:')

    axes[1].plot(t_spikes_2, i_spikes_2 + 1, '.k', markersize=15)
    axes[1].axis([t_final - 10, t_final, 0, N + 1])
    axes[1].set_yticks([1, 2, 3])
    axes[1].set_xlabel('$t$ [units of $T$]')
    axes[1].set_title(r'$\delta=0.55$:')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
