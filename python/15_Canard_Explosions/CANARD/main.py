import numpy as np
import matplotlib.pyplot as plt

a = 5.
tau_n = 60.

t_final = 10000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)

tail_start = 4 * m_steps // 5


def _amplitude(I):
    '''Heun/RK2 integration using plain floats (no array storage) -- just
    the max-min of v over the tail 1/5 of the run, needed by the bisection
    below. Plain-float arithmetic here instead of numpy scalars/arrays is
    what keeps ~100 of these (1e6 steps each) practical to run.'''
    v, n = 0., 0.1
    v_max, v_min = -np.inf, np.inf
    for k in range(m_steps):
        v_inc = v - v ** 3 / 3 - n + I
        n_inc = (a * v - n) / tau_n
        v_tmp = v + dt05 * v_inc
        n_tmp = n + dt05 * n_inc
        v_inc = v_tmp - v_tmp ** 3 / 3 - n_tmp + I
        n_inc = (a * v_tmp - n_tmp) / tau_n
        v = v + dt * v_inc
        n = n + dt * n_inc
        if k >= tail_start:
            if v > v_max:
                v_max = v
            if v < v_min:
                v_min = v
    return v_max - v_min


def simulate(I):
    '''same integration, but returns the full (v, n) trajectories'''
    v = np.zeros(m_steps + 1)
    n = np.zeros(m_steps + 1)
    v[0], n[0] = 0., 0.1
    for k in range(m_steps):
        v_inc = v[k] - v[k] ** 3 / 3 - n[k] + I
        n_inc = (a * v[k] - n[k]) / tau_n
        v_tmp = v[k] + dt05 * v_inc
        n_tmp = n[k] + dt05 * n_inc
        v_inc = v_tmp - v_tmp ** 3 / 3 - n_tmp + I
        n_inc = (a * v_tmp - n_tmp) / tau_n
        v[k + 1] = v[k] + dt * v_inc
        n[k + 1] = n[k] + dt * n_inc
    return v, n


def find_I_for_amplitude(amp_target, i_left=-4.28, i_right=-4.25, tol=1e-10):
    '''bisect on I so the limit-cycle amplitude (max-min of v, tail 1/5 of
    the run) matches amp_target -- the canard explosion is so steep that
    the whole amplitude range 1..3.78 lives in a tiny window of I'''
    while i_right - i_left > tol:
        I = (i_left + i_right) / 2
        amp = _amplitude(I)
        if amp > amp_target:
            i_right = I
        else:
            i_left = I
    I = (i_left + i_right) / 2
    return I, simulate(I)


amp_vec = [1, 2, 3, 3.5, 3.78]
results = []
for amp_target in amp_vec:
    I, (v, n) = find_I_for_amplitude(amp_target)
    results.append((amp_target, I, v, n))

if __name__ == "__main__":

    fig = plt.figure(figsize=(10, 8))
    ax_phase = fig.add_subplot(2, 2, 1)
    ax_t = {1: fig.add_subplot(2, 2, 2), 2: fig.add_subplot(2, 2, 3),
            5: fig.add_subplot(2, 2, 4)}

    A, B, C, D = -3, 2, -6, -2
    L_200 = round(200 / dt)

    for ijk, (amp_target, I, v, n) in enumerate(results, start=1):
        half = v[m_steps // 2:]
        n_half = n[m_steps // 2:]
        if amp_target == 3.5:
            ax_phase.plot(half, n_half, '-b', linewidth=3)
        else:
            ax_phase.plot(half, n_half, '-k', linewidth=1)

        if ijk in ax_t:
            t = np.arange(L_200 + 1) * dt
            ax_t[ijk].plot(t, half[-(L_200 + 1):], '-k', linewidth=2)
            ax_t[ijk].set_xlabel('$t$')
            ax_t[ijk].set_ylabel('$v$')
            ax_t[ijk].set_title(f'$I={I:.7g}$')
            ax_t[ijk].set_xlim(0, 200)
            ax_t[ijk].set_ylim(-3, 2)

    v_line = A + np.arange(101) / 100 * (B - A)
    I_last = results[-1][1]
    ax_phase.plot(v_line, v_line - v_line ** 3 / 3 + I_last, '-g', linewidth=2)
    ax_phase.plot(v_line, a * v_line, '-r', linewidth=2)
    ax_phase.set_xlim(A, B)
    ax_phase.set_ylim(C, D)
    ax_phase.set_box_aspect(1)
    ax_phase.set_xlabel('$v$')
    ax_phase.set_ylabel('$n$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
