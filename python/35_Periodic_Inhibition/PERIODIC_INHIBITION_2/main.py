import numpy as np
import matplotlib.pyplot as plt

# MATLAB's rng('default'); rng(63806) cannot be bit-reproduced by NumPy,
# so we use our own seed and verify statistically/visually instead of
# expecting an exact match to MATLAB's spike times.
rng = np.random.default_rng(63806)

alpha = 5.
Period = 25.
g_bar = 0.1
N = 2000
m = np.mean(np.exp(alpha * np.cos(np.pi * np.arange(N) / N) ** 2) - 1)
tau = 10.
I = 0.1
tau_noise = 3.
sigma_noise = 0.08

t_final = 500.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)
gamma = sigma_noise * np.sqrt(1 - np.exp(-2 * dt / tau_noise))

t = np.arange(m_steps + 1) * dt


def g(t):
    return g_bar * (np.exp(alpha * np.cos(np.pi * t / Period) ** 2) - 1) / m


s_noise = np.zeros(m_steps + 1)
s_noise[0] = sigma_noise * rng.standard_normal()
for k in range(m_steps):
    s_noise[k + 1] = s_noise[k] * np.exp(-dt / tau_noise) + gamma * rng.standard_normal()


def run(use_periodic):
    v = np.zeros(m_steps + 1)
    k_old = 0
    num_spikes = 0
    segments = []
    spike_times = []
    for k in range(1, m_steps + 1):
        g_old = g((k - 1) * dt) if use_periodic else g_bar
        g_mid = g((k - 0.5) * dt) if use_periodic else g_bar
        v_inc = -v[k - 1] / tau + I + s_noise[k - 1] - g_old * v[k - 1]
        v_tmp = v[k - 1] + dt05 * v_inc
        v_inc = -v_tmp / tau + I + (s_noise[k - 1] + s_noise[k]) / 2 - g_mid * v_tmp
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
    fig, axes = plt.subplots(3, 1, figsize=(9, 8))

    axes[0].plot(t, g(t), '-b', linewidth=2)
    axes[0].plot(t, g_bar * np.ones(m_steps + 1), '-r', linewidth=2)
    axes[0].set_ylabel('$g$')
    axes[0].axis([0, t_final, 0, 0.5])

    segs, spikes, freq = run(use_periodic=True)
    for seg_t, seg_v in segs:
        axes[1].plot(seg_t, seg_v, '-b', linewidth=2)
    for ts in spikes:
        axes[1].plot([ts, ts], [0, 5], '-b', linewidth=2)
    axes[1].axis([0, t_final, -1, 6])
    axes[1].set_ylabel('$v$')
    print(f"freq={freq:.2f} Hz")

    segs_bar, spikes_bar, freq_bar = run(use_periodic=False)
    for seg_t, seg_v in segs_bar:
        axes[2].plot(seg_t, seg_v, '-r', linewidth=2)
    for ts in spikes_bar:
        axes[2].plot([ts, ts], [0, 5], '-r', linewidth=2)
    axes[2].axis([0, t_final, -1, 6])
    axes[2].set_ylabel(r'$\overline{v}$')
    axes[2].set_xlabel('$t$ [ms]')
    print(f"freq_bar={freq_bar:.2f} Hz")

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
