import numpy as np


def draw_arrow(ax, xlim, ylim, x, y, v, epsilon=0.1, width=2, color='k'):
    """Draw a V-shaped arrowhead at (x, y) pointing back along -v.

    Port of the book's matlab/*/arrow.m helper (used in 16+ chapters).
    v is normalized by the axis data range first, so the arrowhead looks
    right even when x and y are on very different scales, then rotated
    +-30 degrees for the two limbs.
    """
    a, b = xlim
    c, d = ylim
    v = np.asarray(v, dtype=float)
    u = np.array([v[0] / (b - a), v[1] / (d - c)])
    u = u / np.linalg.norm(u) * epsilon

    rot_right = np.array([[np.cos(np.pi / 6), np.sin(np.pi / 6)],
                           [-np.sin(np.pi / 6), np.cos(np.pi / 6)]])
    rot_left = np.array([[np.cos(np.pi / 6), -np.sin(np.pi / 6)],
                          [np.sin(np.pi / 6), np.cos(np.pi / 6)]])

    for rot in (rot_right, rot_left):
        u_rot = rot @ u
        v_rot = np.array([u_rot[0] * (b - a), u_rot[1] * (d - c)])
        ax.plot([x, x - v_rot[0]], [y, y - v_rot[1]], color=color, linewidth=width)
