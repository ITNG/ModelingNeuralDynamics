import numpy as np
import matplotlib.pyplot as plt

tau_m = 10.
I = 0.12
g = 0.15
tau_I = 7.

dt = 0.01
dt05 = dt / 2


def run(v0, continue_while):
    '''Heun-integrate a LIF neuron receiving a decaying inhibitory
    conductance g*exp(-t/tau_I) starting at t=0, from v(0)=v0 until
    continue_while(v_k) is False. Returns the voltage trace and the
    interpolated threshold-crossing time ("period").'''
    v = [v0]
    k = 0
    while continue_while(v[k]):
        t = k * dt
        v_inc = -v[k] / tau_m + I - g * np.exp(-t / tau_I) * v[k]
        v_tmp = v[k] + dt05 * v_inc
        v_inc = -v_tmp / tau_m + I - g * np.exp(-(t + dt05) / tau_I) * v_tmp
        v.append(v[k] + dt * v_inc)
        k += 1
    v = np.array(v)
    period = ((k - 1) * dt * (v[k] - 1) + k * dt * (1 - v[k - 1])) / (v[k] - v[k - 1])
    return v, k, period


v_0, k_0, period_0 = run(0., lambda vk: vk < 1)
v_1, k_1, period_1 = run(1., lambda vk: vk <= 1)

if __name__ == "__main__":

    plt.figure(figsize=(8, 5))

    t_0 = np.arange(k_0) * dt
    plt.plot(t_0, v_0[:k_0], '-k', linewidth=2)
    plt.plot([(k_0 - 1) * dt, period_0], [v_0[k_0 - 1], 1], '-k', linewidth=2)
    plt.plot([period_0, period_0], [0, 1], ':k', linewidth=2)
    plt.text(period_0 - 0.5, -0.15, '$P_0$', fontsize=16)

    t_1 = np.arange(k_1) * dt
    plt.plot(t_1, v_1[:k_1], '-k', linewidth=2)
    plt.plot([(k_1 - 1) * dt, period_1], [v_1[k_1 - 1], 1], '-k', linewidth=2)
    plt.plot([period_1, period_1], [0, 1], ':k', linewidth=2)
    plt.text(period_1 - 1, -0.15, '$P_1$', fontsize=16)

    plt.xticks(range(0, 31, 15))
    plt.axis([0, 30, 0, 1])
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
