import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from mnd.core import draw_arrow

# The trajectories plotted here, and the locations of the arrows, were
# put in with a lot of trial and error to make the figure look nice.
# Don't try to detect a system, there was none. (matlab source's own
# comment.)

tau_m = 1.
I = 0.3
tau_I = 9.

dt = 0.01
dt05 = dt / 2


def rhs(theta, g_syn):
    theta_inc = -np.cos(theta) / tau_m + (2 * I - g_syn) * (1 + np.cos(theta)) - g_syn * np.sin(theta)
    g_syn_inc = -g_syn / tau_I
    return theta_inc, g_syn_inc


def trajectory(theta0, g_syn0, forward=True):
    theta = [theta0]
    g_syn = [g_syn0]
    k = 0
    sign = 1 if forward else -1
    while -np.pi <= theta[k] <= np.pi:
        theta_inc, g_syn_inc = rhs(theta[k], g_syn[k])
        theta_tmp = theta[k] + sign * dt05 * theta_inc
        g_syn_tmp = g_syn[k] + sign * dt05 * g_syn_inc
        theta_inc, g_syn_inc = rhs(theta_tmp, g_syn_tmp)
        theta.append(theta[k] + sign * dt * theta_inc)
        g_syn.append(g_syn[k] + sign * dt * g_syn_inc)
        k += 1
    return np.array(theta), np.array(g_syn)


M = 6
trajectories = []
for ijk in range(1, M + 1):
    trajectories.append(trajectory(-np.pi, ijk / M, forward=True))
for ijk in range(1, M + 1):
    trajectories.append(trajectory(np.pi, ijk / M, forward=False))
for ijk in range(1, M + 1):
    trajectories.append(trajectory(0., ijk / (M + 1), forward=True))
for ijk in range(1, M + 1):
    trajectories.append(trajectory(0., ijk / (M + 1), forward=False))

extra_trajectories = [
    trajectory(-2., 1., forward=True),
    trajectory(0., 1., forward=True),
    trajectory(np.pi, 0.075, forward=False),
]

arrow_points = [
    (0., 4 / 7),
    (-2., 0.592),
    (-2., 0.145),
    (0.5, 0.203),
    (2.5, 0.18),
    (2.85, 0.52),
]

theta_sep, g_syn_sep = trajectory(-1.23, 1., forward=True)
g_ast = g_syn_sep[-1]
arrow_sep = (-0.75, 0.25)

theta_unstable, g_syn_unstable = trajectory(1.25, 0.25, forward=False)

if __name__ == "__main__":

    fig, ax = plt.subplots(figsize=(7, 7))

    for theta, g_syn in trajectories + extra_trajectories:
        ax.plot(theta, g_syn, '-k', linewidth=2)

    for theta0, g_syn0 in arrow_points:
        v = rhs(theta0, g_syn0)
        draw_arrow(ax, (-np.pi, np.pi), (0, 1), theta0, g_syn0, v, epsilon=0.05, width=2, color='k')

    ax.plot(theta_sep, g_syn_sep, '-b', linewidth=4)
    v = rhs(*arrow_sep)
    draw_arrow(ax, (-np.pi, np.pi), (0, 1), arrow_sep[0], arrow_sep[1], v, epsilon=0.05, width=4, color='b')

    ax.plot(theta_unstable, g_syn_unstable, '-r', linewidth=4)

    ax.plot(np.pi, g_ast, '.k', markersize=15)
    ax.text(np.pi + 0.1, g_ast, r'$g_\ast$', fontsize=20)

    ax.axis([-np.pi, np.pi, 0, 1])
    ax.set_box_aspect(1)
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$g_{\rm syn}$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
