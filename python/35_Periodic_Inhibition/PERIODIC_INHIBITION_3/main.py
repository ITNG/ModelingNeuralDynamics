import numpy as np
import matplotlib.pyplot as plt

tau = 10.
alpha = 1.
Period = 25.
g_bar = 0.1
N = 2000
m = np.mean(np.exp(alpha * np.cos(np.pi * np.arange(N) / N) ** 2) - 1)


def g(t):
    return g_bar * (np.exp(alpha * np.cos(np.pi * t / Period) ** 2) - 1) / m


t_final = 100.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)
t = np.arange(m_steps + 1) * dt


def run(I, use_periodic):
    '''Integrate the LIF neuron; return (t, v) segments split at each
    reset (matlab plots each inter-spike segment separately), the
    spike times (for the vertical reset markers), and the firing
    frequency.'''
    v = np.zeros(m_steps + 1)
    k_old = 0
    num_spikes = 0
    segments = []
    spike_times = []
    for k in range(1, m_steps + 1):
        g_old = g((k - 1) * dt) if use_periodic else g_bar
        g_mid = g((k - 0.5) * dt) if use_periodic else g_bar
        v_inc = -v[k - 1] / tau + I - g_old * v[k - 1]
        v_tmp = v[k - 1] + dt05 * v_inc
        v_inc = -v_tmp / tau + I - g_mid * v_tmp
        v[k] = v[k - 1] + dt * v_inc
        if v[k] > 1:
            segments.append((t[k_old:k + 1], v[k_old:k + 1].copy()))
            k_old = k
            v[k] = 0
            num_spikes += 1
            spike_times.append(k * dt)
    segments.append((t[k_old:], v[k_old:]))
    freq = num_spikes / t_final * 1000
    return segments, spike_times, freq


if __name__ == "__main__":
    fig, axes = plt.subplots(3, 1, figsize=(7, 8))

    axes[0].plot(t, g(t), '-b', linewidth=2)
    axes[0].plot([0, t_final], [g_bar, g_bar], '-r', linewidth=2)
    axes[0].set_title(r'$g$ (blue) and $\overline{g}$ (red)')
    axes[0].axis([0, t_final, 0, g(t).max() * 1.2])

    for ax, I in zip(axes[1:], [0.15, 0.2]):
        segs, spikes, freq = run(I, use_periodic=True)
        for seg_t, seg_v in segs:
            ax.plot(seg_t, seg_v, '-b', linewidth=2)
        for ts in spikes:
            ax.plot([ts, ts], [0, 5], '-b', linewidth=2)
        segs_bar, spikes_bar, freq_bar = run(I, use_periodic=False)
        for seg_t, seg_v in segs_bar:
            ax.plot(seg_t, seg_v, '-r', linewidth=2)
        for ts in spikes_bar:
            ax.plot([ts, ts], [0, 1], '-r', linewidth=2)
        ax.axis([0, t_final, 0, 6])
        ax.set_title(rf'$v$ (blue) and $\overline{{v}}$ (red), $I={I}$')
        print(f"I={I}: freq={freq:.2f} Hz, freq_bar={freq_bar:.2f} Hz")

    axes[-1].set_xlabel('$t$ [ms]')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
