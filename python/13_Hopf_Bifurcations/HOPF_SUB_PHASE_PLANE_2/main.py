import numpy as np
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

dt = 0.001


def spiral(I, r0, theta0, t_final, reverse=False):
    '''dr/dt = I*r + r^3 - r^5 (quintic Hopf normal form radial dynamics).

    reverse=True traces the trajectory backward in time (same idea as
    HOPF_SUB_PHASE_PLANE's spiral()) -- numerically stable way to draw
    the part of the orbit near a fixed point or limit cycle that a
    forward integration would take unboundedly long to leave/reach.
    '''
    sign = -1 if reverse else 1
    m_steps = round(t_final / dt)
    r, theta = np.zeros(m_steps + 1), np.zeros(m_steps + 1)
    r[0], theta[0] = r0, theta0
    for k in range(m_steps):
        r_inc = sign * (I * r[k] + r[k] ** 3 - r[k] ** 5)
        r_tmp = r[k] + dt / 2 * r_inc
        r_inc = sign * (I * r_tmp + r_tmp ** 3 - r_tmp ** 5)
        r[k + 1] = r[k] + dt * r_inc
        theta[k + 1] = theta[k] + sign * dt
    return r * np.cos(theta), r * np.sin(theta)


def trajectory_with_arrow(ax, xlim, ylim, I, r0, theta0, t_final, phi):
    x, y = spiral(I, r0, theta0, t_final)
    ax.plot(x, y, color='k', linewidth=1)
    u = np.array([x[1] - x[0], y[1] - y[0]])
    rot = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
    u = rot @ u
    draw_arrow(ax, xlim, ylim, x[0], y[0], u, epsilon=0.075, width=1)

    x, y = spiral(I, r0, theta0, t_final, reverse=True)
    ax.plot(x, y, color='k', linewidth=1)


if __name__ == "__main__":

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

    # panel 1: I < -1/4 -- only the origin is a fixed point, all orbits decay
    ax = axes[0]
    A, B, C, D = -.5, .5, -.5, .5
    ax.plot(0, 0, '.k', markersize=18)
    trajectory_with_arrow(ax, [A, B], [C, D], I=-0.3, r0=0.2, theta0=0, t_final=5, phi=0.)
    ax.set_xlim(A, B)
    ax.set_ylim(C, D)
    ax.set_box_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title(r'$I<-1/4$')

    # panel 2: -1/4 < I < 0 -- unstable circle r0 inside stable cycle R0
    ax = axes[1]
    I = -0.2
    r0_cycle = np.sqrt(1 / 2 - np.sqrt(1 / 4 + I))
    R0_cycle = np.sqrt(1 / 2 + np.sqrt(1 / 4 + I))
    A, B, C, D = -1.1, 1.1, -1.1, 1.1
    ax.plot(0, 0, '.k', markersize=18)
    th = np.arange(201) / 200 * 2 * np.pi
    ax.plot(r0_cycle * np.cos(th), r0_cycle * np.sin(th), color='k', linewidth=2, linestyle='dashed')
    ax.plot(R0_cycle * np.cos(th), R0_cycle * np.sin(th), color='k', linewidth=3)
    trajectory_with_arrow(ax, [A, B], [C, D], I, r0=0.4, theta0=0, t_final=5, phi=0.2)
    trajectory_with_arrow(ax, [A, B], [C, D], I, r0=0.7, theta0=np.pi / 2, t_final=5, phi=0.2)
    trajectory_with_arrow(ax, [A, B], [C, D], I, r0=1.1, theta0=5 * np.pi / 4, t_final=5, phi=0.)
    ax.set_xlim(A, B)
    ax.set_ylim(C, D)
    ax.set_box_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title(r'$-1/4<I<0$')

    # panel 3: I > 0 -- only the stable cycle R0 remains, origin unstable
    ax = axes[2]
    I = 0.1
    R0_cycle = np.sqrt(1 / 2 + np.sqrt(1 / 4 + I))
    A, B, C, D = -1.3, 1.3, -1.3, 1.3
    ax.plot(0, 0, 'ok', markerfacecolor='none')
    ax.plot(R0_cycle * np.cos(th), R0_cycle * np.sin(th), color='k', linewidth=3)
    trajectory_with_arrow(ax, [A, B], [C, D], I, r0=1.3, theta0=3 * np.pi / 4, t_final=5, phi=-0.1)
    trajectory_with_arrow(ax, [A, B], [C, D], I, r0=0.5, theta0=np.pi, t_final=5, phi=0.)
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
