import numpy as np
import matplotlib.pyplot as plt

a = 5.
tau_n = 60.


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


i_ext_vec = -6 + np.arange(501) / 500 * 12

branches = {'red': [], 'green': [], 'black': [], 'blue': [], 'magenta': []}
real_part_old = -1000.
i_c = None
for ijk, i_ext in enumerate(i_ext_vec):
    v_c = find_fixed_point(i_ext)
    e = np.linalg.eigvals(jacobian(v_c))

    real_part = e[0].real
    if real_part > 0 and real_part_old < 0:
        # linear interpolation for the I where the fixed point loses stability
        i_c = (i_ext_vec[ijk] * (-real_part_old) + i_ext_vec[ijk - 1] * real_part) \
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

t_final = 200.
dt = 0.01
dt05 = dt / 2
m_steps = round(t_final / dt)
half = m_steps // 2

i_scan = np.arange(-5, 5 + 1e-9, 0.05)
cycle_pts = []
for i in i_scan:
    v, n = 3., 3.
    vmax, vmin = -np.inf, np.inf
    for k in range(m_steps):
        v_inc = v - v ** 3 / 3 - n + i
        n_inc = (a * v - n) / tau_n
        v_tmp = v + dt05 * v_inc
        n_tmp = n + dt05 * n_inc
        v_inc = v_tmp - v_tmp ** 3 / 3 - n_tmp + i
        n_inc = (a * v_tmp - n_tmp) / tau_n
        v = v + dt * v_inc
        n = n + dt * n_inc
        if k >= half:
            if v > vmax:
                vmax = v
            if v < vmin:
                vmin = v
    if vmax - vmin > 0.8:
        cycle_pts.append((i, vmax, vmin))

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

    for i, vmax, vmin in cycle_pts:
        plt.plot(i, vmax, '.k')
        plt.plot(i, vmin, '.k')

    plt.xlabel('$I$')
    plt.ylabel(r'$v_\ast$')
    plt.xlim(i_ext_vec[0], i_ext_vec[-1])
    plt.ylim(-5.7, 4.7)
    plt.gca().set_box_aspect(1)
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
