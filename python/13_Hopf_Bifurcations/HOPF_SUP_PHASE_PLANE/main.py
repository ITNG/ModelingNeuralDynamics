import numpy as np
import matplotlib.pyplot as plt

from mnd.core import draw_arrow


def spiral(I, r0, theta0, t_final, dt, reverse=False):
    '''dr/dt = I*r - r^3 (supercritical Hopf normal form radial dynamics).

    reverse=True traces the trajectory backward in time (same idea as
    HOPF_SUB_PHASE_PLANE's spiral()).
    '''
    sign = -1 if reverse else 1
    m_steps = round(t_final / dt)
    r, theta = np.zeros(m_steps + 1), np.zeros(m_steps + 1)
    r[0], theta[0] = r0, theta0
    for k in range(m_steps):
        r_inc = sign * (I * r[k] - r[k] ** 3)
        r_tmp = r[k] + dt / 2 * r_inc
        r_inc = sign * (I * r_tmp - r_tmp ** 3)
        r[k + 1] = r[k] + dt * r_inc
        theta[k + 1] = theta0 + sign * (k + 1) * dt
    return r * np.cos(theta), r * np.sin(theta)


if __name__ == "__main__":

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # panel 1: I < 0 -- only the origin is a fixed point, all orbits decay
    ax = axes[0]
    A, B, C, D = -.5, .5, -.5, .5
    I = -0.02
    ax.plot(0, 0, '.k', markersize=18)
    x, y = spiral(I, 0.2, 0., 5, dt=0.01)
    ax.plot(x, y, color='k', linewidth=1)
    u = np.array([x[1] - x[0], y[1] - y[0]])
    phi = 0.1
    rot = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
    u = rot @ u
    draw_arrow(ax, [A, B], [C, D], x[0], y[0], u, epsilon=0.075, width=1)
    x, y = spiral(I, 0.2, 0., 5, dt=0.01, reverse=True)
    ax.plot(x, y, color='k', linewidth=1)
    ax.set_xlim(A, B)
    ax.set_ylim(C, D)
    ax.set_box_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title('$I<0$')

    # panel 2: I > 0 -- a stable limit cycle at r=sqrt(I)
    ax = axes[1]
    A, B, C, D = -2, 2, -2, 2
    I = 0.5
    ax.plot(0, 0, 'ok', markerfacecolor='none')
    th = np.arange(101) / 100 * 2 * np.pi
    ax.plot(np.sqrt(I) * np.cos(th), np.sqrt(I) * np.sin(th), color='k', linewidth=4)

    x, y = spiral(I, 0.4, 0., 3, dt=0.001)
    ax.plot(x, y, color='k', linewidth=1)
    u = np.array([x[999] - x[998], y[999] - y[998]])
    phi = 0.3
    rot = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
    u = rot @ u
    draw_arrow(ax, [A, B], [C, D], x[999], y[999], u, epsilon=0.075, width=1)

    x, y = spiral(I, 0.4, 0., 1, dt=0.01, reverse=True)
    ax.plot(x, y, color='k', linewidth=1)

    x, y = spiral(I, 1.2, -np.pi / 2, 3, dt=0.0001)
    ax.plot(x, y, color='k', linewidth=1)
    u = np.array([x[1] - x[0], y[1] - y[0]])
    draw_arrow(ax, [A, B], [C, D], x[0], y[0], u, epsilon=0.075, width=1)

    x, y = spiral(I, 1.2, -np.pi / 2, 3, dt=0.0001, reverse=True)
    ax.plot(x, y, color='k', linewidth=1)

    ax.set_xlim(A, B)
    ax.set_ylim(C, D)
    ax.set_box_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title('$I>0$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
