import numpy as np
import matplotlib.pyplot as plt

a = 5.
tau_n = 60.

i_c = -4.291556094395530
i_ext_vec = 1.01 * i_c + np.arange(501) / 500 * (0.98 * i_c - 1.01 * i_c)


def f(v, i_ext):
    '''zeros of this function are fixed points of the FN system'''
    return v - v ** 3 / 3 - a * v + i_ext


def find_fixed_point(i_ext):
    v_left, v_right = -1., 1.
    while f(v_left, i_ext) < 0:
        v_left -= 1
    while f(v_right, i_ext) > 0:
        v_right += 1
    while v_right - v_left > 1e-10:
        v_c = (v_left + v_right) / 2
        if f(v_c, i_ext) >= 0:
            v_left = v_c
        else:
            v_right = v_c
    return (v_left + v_right) / 2


def jacobian(v_c):
    return np.array([[1 - v_c ** 2, -a], [a / tau_n, -1 / tau_n]])


t_final = 10000.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)
tail_start = 4 * m_steps // 5 - 1  # matlab's v(4*m_steps/5 : m_steps+1) is 1-indexed


def cycle_amplitude(v_c, n_c, i_ext):
    '''plain-float Heun integration (no array storage) -- just the max/min
    of v over the tail 1/5 of a 10s run, starting just outside the
    unstable spiral, needed for the limit-cycle envelope below'''
    v, n = v_c, n_c * 1.05
    vmax, vmin = -np.inf, np.inf
    for k in range(m_steps):
        v_inc = v - v ** 3 / 3 - n + i_ext
        n_inc = (a * v - n) / tau_n
        v_tmp = v + dt05 * v_inc
        n_tmp = n + dt05 * n_inc
        v_inc = v_tmp - v_tmp ** 3 / 3 - n_tmp + i_ext
        n_inc = (a * v_tmp - n_tmp) / tau_n
        v = v + dt * v_inc
        n = n + dt * n_inc
        if k >= tail_start:
            if v > vmax:
                vmax = v
            if v < vmin:
                vmin = v
    return vmax, vmin


branches = {'red': [], 'green': [], 'black': [], 'blue': [], 'magenta': []}
cycle = []  # (i_ext, max_cycle, min_cycle)
real_part_old = -1000.
i_c_found = None

for ijk, i_ext in enumerate(i_ext_vec):
    v_c = find_fixed_point(i_ext)
    n_c = a * v_c
    e = np.linalg.eigvals(jacobian(v_c))

    real_part = e[0].real
    if real_part > 0 and real_part_old < 0:
        i_c_found = (i_ext_vec[ijk] * (-real_part_old) + i_ext_vec[ijk - 1] * real_part) \
            / (real_part - real_part_old)
    real_part_old = real_part

    if abs(e[0].imag) > 1e-4:
        if e[0].real < 0:
            branches['red'].append((i_ext, v_c))
        if e[0].real > 0:
            branches['green'].append((i_ext, v_c))
    else:
        if e[0].real < 0 and e[1].real < 0:
            branches['black'].append((i_ext, v_c))
        if e[0].real > 0 and e[1].real > 0:
            branches['blue'].append((i_ext, v_c))
        if e[0].real > 0 and e[1].real < 0:
            branches['magenta'].append((i_ext, v_c))

    if abs(e[0].imag) > 1e-4 and e[0].real > 0:
        vmax, vmin = cycle_amplitude(v_c, n_c, i_ext)
        cycle.append((i_ext, vmax, vmin))

if __name__ == "__main__":

    plt.figure(figsize=(7, 7))

    if branches['blue']:
        i_pts, v_pts = zip(*branches['blue'])
        plt.plot(i_pts, v_pts, '-.b', linewidth=4)
    if branches['green']:
        i_pts, v_pts = zip(*branches['green'])
        plt.plot(i_pts, v_pts, '--g', linewidth=4)
    if branches['red']:
        i_pts, v_pts = zip(*branches['red'])
        plt.plot(i_pts, v_pts, '.r', markersize=10)
    if branches['black']:
        i_pts, v_pts = zip(*branches['black'])
        plt.plot(i_pts, v_pts, '-k', linewidth=4)
    if branches['magenta']:
        i_pts, v_pts = zip(*branches['magenta'])
        plt.plot(i_pts, v_pts, '--m', linewidth=4)
    if cycle:
        i_pts, max_pts, min_pts = zip(*cycle)
        plt.plot(i_pts, max_pts, ':k', linewidth=2)
        plt.plot(i_pts, min_pts, ':k', linewidth=2)

    plt.xlabel('$I$')
    plt.ylabel(r'$v_\ast$')
    plt.xticks([-4.3, -4.26, -4.22])
    plt.xlim(i_ext_vec[0], i_ext_vec[-1])
    plt.ylim(-2.5, 2)
    plt.gca().set_box_aspect(1)
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
