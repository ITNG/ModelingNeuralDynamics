import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

tau_m = 2.
I = 0.14
T = np.pi * tau_m / np.sqrt(tau_m * I - 0.25)
t_final = 10000.
dt = 0.01
m_steps = round(t_final / dt)
dt05 = dt / 2


def theta_inc(theta):
    return -np.cos(theta) / tau_m + 2 * I * (1 + np.cos(theta))


def simulate(delta, epsilon):
    '''two theta neurons A, B coupled by delayed pulses: whenever B (A)
    spikes, its effect on A (B) arrives delta*T time units later, and
    jumps A's (B's) theta by the phase-response map
    theta -> 2*atan(tan(theta/2) + 2*dv).

    Faithfully reproduces the matlab source's bug in the "no delayed
    input this step" branch for theta_B, which advances theta_B using
    theta_A_inc (leftover from theta_A's update earlier in the same
    step) instead of computing/using its own increment.'''
    dv = np.sqrt(tau_m * I - 0.25) * epsilon

    theta_A = 0.
    theta_B = np.pi / 6
    t_spikes_A = []
    t_spikes_B = []

    # monotonically-advancing pointers into the other neuron's spike
    # list: since delta*T is a fixed offset and spike times only
    # increase, the next not-yet-delivered spike is found in O(1)
    # amortized instead of re-scanning the whole list every step.
    ptr_A = 0  # next unconsumed entry in t_spikes_B (deliveries to A)
    ptr_B = 0  # next unconsumed entry in t_spikes_A (deliveries to B)

    for k in range(1, m_steps + 1):
        t_lo, t_hi = (k - 1) * dt, k * dt

        theta_A_inc = None
        if ptr_A < len(t_spikes_B) and t_lo < t_spikes_B[ptr_A] + delta * T <= t_hi:
            t_0 = t_spikes_B[ptr_A] + delta * T
            ptr_A += 1
            dt_1 = t_0 - t_lo
            dt_105 = dt_1 / 2
            theta_A_inc = theta_inc(theta_A)
            theta_A_tmp = theta_A + dt_105 * theta_A_inc
            theta_A_inc = theta_inc(theta_A_tmp)
            t_A = theta_A + dt_1 * theta_A_inc
            t_A = 2 * np.arctan(np.tan(t_A / 2) + 2 * dv)
            dt_2 = t_hi - t_0
            dt_205 = dt_2 / 2
            theta_A_inc = theta_inc(t_A)
            theta_A_tmp = t_A + dt_205 * theta_A_inc
            theta_A_inc = theta_inc(theta_A_tmp)
            theta_A_next = t_A + dt_2 * theta_A_inc
        else:
            theta_A_inc = theta_inc(theta_A)
            theta_A_tmp = theta_A + dt05 * theta_A_inc
            theta_A_inc = theta_inc(theta_A_tmp)
            theta_A_next = theta_A + dt * theta_A_inc

        if ptr_B < len(t_spikes_A) and t_lo < t_spikes_A[ptr_B] + delta * T <= t_hi:
            t_0 = t_spikes_A[ptr_B] + delta * T
            ptr_B += 1
            dt_1 = t_0 - t_lo
            dt_105 = dt_1 / 2
            theta_B_inc = theta_inc(theta_B)
            theta_B_tmp = theta_B + dt_105 * theta_B_inc
            theta_B_inc = theta_inc(theta_B_tmp)
            t_B = theta_B + dt_1 * theta_B_inc
            t_B = 2 * np.arctan(np.tan(t_B / 2) + 2 * dv)
            dt_2 = t_hi - t_0
            dt_205 = dt_2 / 2
            theta_B_inc = theta_inc(t_B)
            theta_B_tmp = t_B + dt_205 * theta_B_inc
            theta_B_inc = theta_inc(theta_B_tmp)
            theta_B_next = t_B + dt_2 * theta_B_inc
        else:
            # faithful port of the matlab bug: uses theta_A_inc (left
            # over from the A update above), not a freshly computed
            # theta_B_inc.
            theta_B_next = theta_B + dt * theta_A_inc

        if theta_A_next > np.pi:
            tt = (t_lo * (theta_A_next - np.pi) + t_hi * (np.pi - theta_A)) / (theta_A_next - theta_A)
            t_spikes_A.append(tt)
            theta_A_next -= 2 * np.pi

        if theta_B_next > np.pi:
            tt = (t_lo * (theta_B_next - np.pi) + t_hi * (np.pi - theta_B)) / (theta_B_next - theta_B)
            t_spikes_B.append(tt)
            theta_B_next -= 2 * np.pi

        theta_A, theta_B = theta_A_next, theta_B_next

    return np.array(t_spikes_A), np.array(t_spikes_B)


def sync_measure(t_spikes_A, t_spikes_B):
    d = np.array([np.min(np.abs(ta - t_spikes_B)) for ta in t_spikes_A])
    return d[len(t_spikes_A) - 2]


epsilon_vec = np.arange(1, 1001) / 1000 * 5
delta_vec = 1 / 2 - 1 / np.pi * np.arctan(epsilon_vec / 2)

N = 10
grid_points = []  # (epsilon, delta, is_synced)
for i in range(1, N):
    for j in range(1, N):
        delta = j / N
        epsilon = i / N * 5
        t_spikes_A, t_spikes_B = simulate(delta, epsilon)
        synced = sync_measure(t_spikes_A, t_spikes_B) < 1e-2
        grid_points.append((epsilon, delta, synced))

if __name__ == "__main__":

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(epsilon_vec, delta_vec, '-k', linewidth=5)
    ax.set_xlabel(r'$\epsilon$')
    ax.set_ylabel(r'$\delta$')
    ax.axis([0, 5, 0, 1])
    ax.set_box_aspect(1)

    poly_x = np.concatenate(([0.], epsilon_vec, [5.], epsilon_vec[::-1], [0.]))
    poly_y = np.concatenate(([0.5], delta_vec, [1.], np.ones(1000), [1.]))
    ax.add_patch(Polygon(np.column_stack([poly_x, poly_y]), closed=True,
                          facecolor='y', alpha=0.2))

    for epsilon, delta, synced in grid_points:
        if synced:
            ax.plot(epsilon, delta, '.r', markersize=15)
        else:
            ax.plot(epsilon, delta, '.b', markersize=15)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
