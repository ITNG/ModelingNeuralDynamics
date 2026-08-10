import numpy as np
import matplotlib.pyplot as plt

tau_m = 10.
I = 0.12
T = tau_m * np.log(tau_m * I / (tau_m * I - 1))  # firing period without inhibition

N = 10
dt = 0.001
dt05 = dt / 2


def simulate_one(v0, g_syn, tau_I):
    v = [v0]
    k = 0
    while v[k] < 1:
        t = k * dt
        v_inc = -v[k] / tau_m + I - g_syn * np.exp(-t / tau_I) * v[k]
        v_tmp = v[k] + dt05 * v_inc
        v_inc = -v_tmp / tau_m + I - g_syn * np.exp(-(t + dt05) / tau_I) * v_tmp
        v.append(v[k] + dt * v_inc)
        k += 1
    v = np.array(v)
    period = ((k - 1) * dt * (v[k] - 1) + k * dt * (1 - v[k - 1])) / (v[k] - v[k - 1])
    return v, k, period


def simulate_panel(g_syn, tau_I):
    traces = []
    for i in range(1, N + 1):
        v0 = (1 - np.exp(-(N - i) * T / (N * tau_m))) * tau_m * I
        v, k, period = simulate_one(v0, g_syn, tau_I)
        traces.append((v, k, period))
    return traces


traces_none = simulate_panel(g_syn=0., tau_I=1.)  # g_syn=0 -> no inhibition
traces_weak = simulate_panel(g_syn=0.15, tau_I=9.)
traces_strong = simulate_panel(g_syn=2., tau_I=1.)

if __name__ == "__main__":

    fig, axes = plt.subplots(3, 1, figsize=(8, 9))

    for ax, traces in zip(axes, (traces_none, traces_weak, traces_strong)):
        for v, k, period in traces:
            t = np.arange(k) * dt
            ax.plot(t, v[:k], '-k', linewidth=1)
            ax.plot([(k - 1) * dt, period], [v[k - 1], 1], '-k', linewidth=1)
            ax.plot([period, period], [0, 1], ':k', linewidth=1)
        ax.axis([0, 30, 0, 1])

    axes[-1].set_xlabel('$t$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
