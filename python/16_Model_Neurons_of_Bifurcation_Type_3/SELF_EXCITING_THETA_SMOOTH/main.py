import numpy as np
import matplotlib.pyplot as plt

I = -0.05
w_max = 0.2
tau_w = 20.

t_final = 50.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

theta = np.zeros(m_steps + 1)
w = np.zeros(m_steps + 1)
theta[0] = np.pi / 2

for k in range(m_steps):
    theta_inc = 1 - np.cos(theta[k]) + (I + w[k]) * (1 + np.cos(theta[k]))
    w_inc = -w[k] / tau_w + 10 * np.exp(-5 * (1 + np.cos(theta[k]))) * (w_max - w[k])
    theta_tmp = theta[k] + dt05 * theta_inc
    w_tmp = w[k] + dt05 * w_inc
    theta_inc = 1 - np.cos(theta_tmp) + (I + w_tmp) * (1 + np.cos(theta_tmp))
    w_inc = -w_tmp / tau_w + 10 * np.exp(-5 * (1 + np.cos(theta_tmp))) * (w_max - w_tmp)
    theta[k + 1] = theta[k] + dt * theta_inc
    w[k + 1] = w[k] + dt * w_inc

if __name__ == "__main__":

    t = np.arange(m_steps + 1) * dt

    fig, ax = plt.subplots(2, figsize=(7, 6), sharex=True)
    ax[0].plot(t, 1 - np.cos(theta), '-k', linewidth=2)
    ax[0].set_xlim(0, t_final)
    ax[0].set_ylim(0, 2)
    ax[0].set_ylabel(r'$1-\cos\theta$')

    ax[1].plot(t, w, '-k', linewidth=2)
    if w_max > 0:
        ax[1].set_ylim(0, w_max)
    ax[1].set_xlabel('$t$ [ms]')
    ax[1].set_ylabel('$z$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
