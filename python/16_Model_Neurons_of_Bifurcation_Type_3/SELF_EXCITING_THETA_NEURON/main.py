import numpy as np
import matplotlib.pyplot as plt

I = -0.05
w_max = 0.20
tau_w = 20.

t_final = 50.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

theta = np.zeros(m_steps + 1)
w = np.zeros(m_steps + 1)
theta[0] = np.pi / 2
k_spikes = []

for k in range(m_steps):
    theta_inc = 1 - np.cos(theta[k]) + (I + w[k]) * (1 + np.cos(theta[k]))
    w_inc = -w[k] / tau_w
    theta_tmp = theta[k] + dt05 * theta_inc
    w_tmp = w[k] + dt05 * w_inc
    theta_inc = 1 - np.cos(theta_tmp) + (I + w_tmp) * (1 + np.cos(theta_tmp))
    w_inc = -w_tmp / tau_w
    theta[k + 1] = theta[k] + dt * theta_inc
    w[k + 1] = w[k] + dt * w_inc
    if theta[k + 1] > np.pi:
        theta[k + 1] = -np.pi
        w[k + 1] = w_max
        k_spikes.append(k)

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    fig, ax = plt.subplots(2, figsize=(7, 6), sharex=True)
    ax[0].plot(t, 1 - np.cos(theta), '-k', linewidth=2)
    ax[0].set_ylim(0, 2)
    ax[0].set_ylabel(r'$1-\cos\theta$')

    ax[1].plot(t[:k_spikes[0] + 1], w[:k_spikes[0] + 1], '-k', linewidth=2)
    t_mid = (t[k_spikes[0]] + t[k_spikes[0] + 1]) / 2
    ax[1].plot([t_mid, t_mid], [w[k_spikes[0]], w_max], ':k', linewidth=2)
    for i in range(len(k_spikes) - 1):
        ax[1].plot(t[k_spikes[i] + 1:k_spikes[i + 1] + 1],
                    w[k_spikes[i] + 1:k_spikes[i + 1] + 1], '-k', linewidth=2)
        t_mid = (t[k_spikes[i + 1]] + t[k_spikes[i + 1] + 1]) / 2
        ax[1].plot([t_mid, t_mid], [w[k_spikes[i + 1]], w_max], ':k', linewidth=2)
    ax[1].plot(t[k_spikes[-1] + 1:], w[k_spikes[-1] + 1:], '-k', linewidth=2)

    ax[1].set_xlim(0, t_final)
    if w_max > 0:
        ax[1].set_ylim(0, w_max)
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$z$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
