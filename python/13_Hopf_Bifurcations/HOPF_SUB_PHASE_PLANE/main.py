import numpy as np
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

A, B, C, D = -.5, .5, -.5, .5
dt = 0.01


def spiral(I, r0, theta0, t_final, reverse=False):
    '''dr/dt = I*r + r^3 (subcritical Hopf normal form radial dynamics).

    reverse=True integrates -f(r) forward with theta decreasing instead --
    traces the same trajectory type (unstable circle / manifold) from the
    other temporal direction, which is numerically stable to compute even
    though the true forward trajectory would take unboundedly long to
    leave the neighborhood of a fixed point.
    '''
    sign = -1 if reverse else 1
    m_steps = round(t_final / dt)
    r, theta = np.zeros(m_steps + 1), np.zeros(m_steps + 1)
    r[0], theta[0] = r0, theta0
    for k in range(m_steps):
        r_inc = sign * (I * r[k] + r[k] ** 3)
        r_tmp = r[k] + dt / 2 * r_inc
        r_inc = sign * (I * r_tmp + r_tmp ** 3)
        r[k + 1] = r[k] + dt * r_inc
        theta[k + 1] = theta[k] + sign * dt
    return r * np.cos(theta), r * np.sin(theta)


def panel(ax, I, r0_pairs, title):
    if I < 0:
        ax.plot(0, 0, '.k', markersize=18)
    else:
        ax.plot(0, 0, 'ok', markerfacecolor='none')

    if I < 0:
        th = np.arange(101) / 100 * 2 * np.pi
        ax.plot(np.sqrt(-I) * np.cos(th), np.sqrt(-I) * np.sin(th),
                color='k', linewidth=1, linestyle='dashed')

    for r0, t_fwd, t_rev, phi in r0_pairs:
        x, y = spiral(I, r0, 0., t_fwd)
        ax.plot(x, y, color='k', linewidth=1)
        u = np.array([x[1] - x[0], y[1] - y[0]])
        rot = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
        u = rot @ u
        draw_arrow(ax, [A, B], [C, D], x[0], y[0], u, epsilon=0.075, width=1)

        x, y = spiral(I, r0, 0., t_rev, reverse=True)
        ax.plot(x, y, color='k', linewidth=1)

    ax.set_xlim(A, B)
    ax.set_ylim(C, D)
    ax.set_box_aspect(1)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_title(title)


if __name__ == "__main__":

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    panel(axes[0], I=-0.1,
          r0_pairs=[(0.1, 5, 10, 0.3), (0.4, 3, 3, 0.1)],
          title='$I<0$')
    axes[0].set_xticks([])
    axes[0].set_yticks([])

    panel(axes[1], I=0.1,
          r0_pairs=[(0.3, 5, 10, 0.1)],
          title='$I>0$')
    axes[1].set_xticks([])
    axes[1].set_yticks([])

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
