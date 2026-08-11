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


def alpha_h(v):
    """Traub-Miles HH-type gating rate."""
    return 0.128 * np.exp(-(v + 50.0) / 18.0)


def alpha_m(v):
    return 0.32 * (v + 54) / (1.0 - np.exp(-(v + 54.0) / 4.0))


def alpha_n(v):
    return 0.032 * (v + 52) / (1.0 - np.exp(-(v + 52.0) / 5.0))


def beta_h(v):
    return 4.0 / (1.0 + np.exp(-(v + 27.0) / 5.0))


def beta_m(v):
    return 0.28 * (v + 27.0) / (np.exp((v + 27.0) / 5.0) - 1.0)


def beta_n(v):
    return 0.5 * np.exp(-(v + 57.0) / 40.0)


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))
