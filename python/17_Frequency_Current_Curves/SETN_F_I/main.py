import numpy as np
import matplotlib.pyplot as plt

tau_z = 100.
z_max = 0.05

I_L, I_R = -0.06, 0.02
i_ext_vec = I_L + np.arange(31) / 30 * (I_R - I_L)

dt = 0.01
dt05 = dt / 2
max_num_spikes = 3
t_max = 1000.


def run(i_ext, theta0):
    '''Heun/RK2 integration (plain floats) of the self-exciting theta
    neuron from a fixed initial condition, for up to t_max ms or until
    max_num_spikes spikes are found. Returns the frequency from the last
    two spike times, or 0 if fewer than max_num_spikes spikes occurred.'''
    theta, z = theta0, 0.
    num_spikes = 0
    t_spikes = []
    k = 0
    t = 0.
    while num_spikes < max_num_spikes and t < t_max:
        theta_inc = 1 - np.cos(theta) + (i_ext + z) * (1 + np.cos(theta))
        z_inc = -z / tau_z + 10 * np.exp(-5 * (1 + np.cos(theta))) * (z_max - z)
        theta_tmp = theta + dt05 * theta_inc
        z_tmp = z + dt05 * z_inc
        theta_inc = 1 - np.cos(theta_tmp) + (i_ext + z_tmp) * (1 + np.cos(theta_tmp))
        z_inc = -z_tmp / tau_z + 10 * np.exp(-5 * (1 + np.cos(theta_tmp))) * (z_max - z_tmp)
        theta_next = theta + dt * theta_inc
        z = z + dt * z_inc
        k += 1
        if theta_next > np.pi:
            num_spikes += 1
            t_spike = ((k - 1) * dt * (theta_next - np.pi) + k * dt * (np.pi - theta)) \
                / (theta_next - theta)
            t_spikes.append(t_spike)
            theta_next -= 2 * np.pi
        theta = theta_next
        t += dt

    if num_spikes == max_num_spikes:
        return 1000 / (t_spikes[-1] - t_spikes[-2])
    return 0.


f_low = np.array([run(i_ext, 0.) for i_ext in i_ext_vec])
f_high = np.array([run(i_ext, 9 / 10 * np.pi) for i_ext in i_ext_vec])

if __name__ == "__main__":

    plt.figure(figsize=(7, 3.5))
    plt.plot(i_ext_vec, f_low, '.k', markersize=15)
    plt.plot(i_ext_vec, f_high, 'ok', markersize=10, markerfacecolor='none', linewidth=1)
    plt.xlim(i_ext_vec.min(), i_ext_vec.max())
    plt.ylim(0, 100)
    plt.xlabel('$I$')
    plt.ylabel('$f$')
    plt.xticks(np.arange(-0.06, 0.021, 0.02))
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
