import numpy as np
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

b = 1.0
# domain trajectories start from (some start off-screen and flow inward)
x_min, x_max = 0.2, 0.2 + 3
y_min, y_max = -0.5, 2.5
# the actually-plotted, clipped view
frame_xlim, frame_ylim = [0.2, 2.2], [-0.5, 1.5]

t_final = 200.
dt = 0.01
rad = 0.1
arrow_eps = 0.125 * np.sqrt(5 / 4)


def derivative(x, y, a):
    dx = -a * x + y
    dy = x ** 2 / (1 + x ** 2) - b * y
    return dx, dy


def simulate(x0, y0, a):
    m_steps = round(t_final / dt)
    x, y = np.zeros(m_steps + 1), np.zeros(m_steps + 1)
    x[0], y[0] = x0, y0
    for k in range(m_steps):
        dx, dy = derivative(x[k], y[k], a)
        x_tmp, y_tmp = x[k] + dt / 2 * dx, y[k] + dt / 2 * dy
        dx, dy = derivative(x_tmp, y_tmp, a)
        x[k + 1] = x[k] + dt * dx
        y[k + 1] = y[k] + dt * dy
    return x, y


def trajectory(ax, x0, y0, a, arrow_y0s=()):
    x, y = simulate(x0, y0, a)
    for y0_cross in arrow_y0s:
        for k in np.where((y[1:] - y0_cross) * (y[:-1] - y0_cross) < 0)[0]:
            vec = [x[k + 1] - x[k], y[k + 1] - y[k]]
            draw_arrow(ax, frame_xlim, frame_ylim, x[k + 1], y[k + 1], vec,
                       epsilon=arrow_eps, width=1)
    ax.plot(x, y, color='k', linewidth=1)


def fill_circle(ax, cx, cy, color, theta_start=0., theta_span=2 * np.pi):
    theta = np.linspace(theta_start, theta_start + theta_span, 101)
    ax.fill(cx + rad * np.cos(theta), cy + rad * np.sin(theta), color=color,
            edgecolor='k' if color == 'w' else None, linewidth=1)


def frame(ax):
    ax.plot([frame_xlim[0], frame_xlim[1], frame_xlim[1], frame_xlim[0], frame_xlim[0]],
            [frame_ylim[0], frame_ylim[0], frame_ylim[1], frame_ylim[1], frame_ylim[0]],
            color='k', linewidth=2)
    ax.set_xlim(*frame_xlim)
    ax.set_ylim(*frame_ylim)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])


def panel(ax, a, left, bottom, right, top, fixed_points):
    if fixed_points:
        disc = 1 / (4 * a ** 2 * b ** 2) - 1
        x_plus = 1 / (2 * a * b) + np.sqrt(disc)
        x_minus = 1 / (2 * a * b) - np.sqrt(disc)
        fill_circle(ax, x_plus, a * x_plus, 'k')

    for edges in (left, bottom, right, top):
        for x0, y0, arrows in edges:
            trajectory(ax, x0, y0, a, arrows)

    frame(ax)

    if fixed_points:
        if abs(a * b - 0.5) < 1e-9:
            # exactly at bifurcation: x_plus == x_minus, draw the white
            # half over the black circle already there -> half-stable mark
            fill_circle(ax, x_minus, a * x_minus, 'w', theta_start=0.73 * np.pi, theta_span=np.pi)
        else:
            fill_circle(ax, x_minus, a * x_minus, 'w')


def edge_left(offset_ijk=4, offset=0.05, arrows_by_ijk=None):
    arrows_by_ijk = arrows_by_ijk or {}
    out = []
    for ijk in range(11):
        y0 = y_min + ijk * (y_max - y_min) / 10
        if ijk == offset_ijk:
            y0 += offset
        out.append((x_min, y0, arrows_by_ijk.get(ijk, ())))
    return out


def edge_bottom(offset_ijk=4, offset=-0.05, arrows_by_ijk=None):
    arrows_by_ijk = arrows_by_ijk or {}
    out = []
    for ijk in range(11):
        x0 = x_min + ijk * (x_max - x_min) / 10
        if ijk == offset_ijk:
            x0 += offset
        out.append((x0, y_min, arrows_by_ijk.get(ijk, ())))
    return out


def edge_right(arrows_by_ijk=None):
    arrows_by_ijk = arrows_by_ijk or {}
    return [(x_max, y_min + ijk * (y_max - y_min) / 10, arrows_by_ijk.get(ijk, ()))
            for ijk in range(4)]


def edge_top(arrows_by_ijk=None):
    arrows_by_ijk = arrows_by_ijk or {}
    return [(x_min + ijk * (x_max - x_min) / 10, y_max, arrows_by_ijk.get(ijk, ()))
            for ijk in range(4)]


if __name__ == "__main__":

    fig, axes = plt.subplots(1, 3, figsize=(13, 5))

    # panel 1: a=0.45 -- two fixed points (stable node + saddle)
    panel(axes[0], a=0.45,
          left=edge_left(arrows_by_ijk={4: [0.5, 0.08]}),
          bottom=edge_bottom(arrows_by_ijk={4: [0, 0.58], 7: [0.0]}),
          right=edge_right(arrows_by_ijk={1: [0.65]}),
          top=edge_top(arrows_by_ijk={0: [1.0]}),
          fixed_points=True)

    # panel 2: a=0.5 -- exactly at the bifurcation, one half-stable point
    panel(axes[1], a=0.5,
          left=edge_left(arrows_by_ijk={5: [0.3], 6: [0.85]}),
          bottom=edge_bottom(arrows_by_ijk={6: [0], 10: [0.35]}),
          right=edge_right(),
          top=edge_top(arrows_by_ijk={0: [1.0, 0.61]}),
          fixed_points=True)

    # panel 3: a=0.55 -- past the bifurcation, no fixed points left
    panel(axes[2], a=0.55,
          left=edge_left(arrows_by_ijk={6: [0.75], 7: [0.375], 8: [1.0]}),
          bottom=edge_bottom(arrows_by_ijk={5: [0], 8: [0]}),
          right=edge_right(),
          top=edge_top(),
          fixed_points=False)

    for ax in axes:
        ax.set_xlabel('$x$')
    axes[0].set_ylabel('$y$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
