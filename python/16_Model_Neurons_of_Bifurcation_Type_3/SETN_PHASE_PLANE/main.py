import numpy as np
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

w_max = 0.2
tau_w = 20.
dt = 0.0001
dt05 = dt / 2

A, B, C, D = -np.pi, np.pi, 0., 2 * w_max


def simulate(I, w0):
    '''Heun/RK2 integration (plain floats/lists) of the self-exciting theta
    neuron until theta wraps past pi or w decays to ~0'''
    theta = [-np.pi]
    w = [w0]
    while theta[-1] < np.pi and w[-1] > 1e-6:
        th, wk = theta[-1], w[-1]
        theta_inc = 1 - np.cos(th) + (I + wk) * (1 + np.cos(th))
        w_inc = -wk / tau_w
        theta_tmp = th + dt05 * theta_inc
        w_tmp = wk + dt05 * w_inc
        theta_inc = 1 - np.cos(theta_tmp) + (I + w_tmp) * (1 + np.cos(theta_tmp))
        w_inc = -w_tmp / tau_w
        theta.append(th + dt * theta_inc)
        w.append(wk + dt * w_inc)
    return np.array(theta), np.array(w)


def theta_thresholds(I):
    if I <= 0:
        theta_plus = np.arccos((1 + I) / (1 - I))
    else:
        theta_plus = 10.
    return theta_plus, -theta_plus


def _arrow(ax, theta, w, ind, width):
    x, y = theta[ind], w[ind]
    vec = np.array([theta[ind + 1] - theta[ind], w[ind + 1] - w[ind]])
    draw_arrow(ax, (A, B), (C, D), x, y, vec, epsilon=0.06, width=width, color='k')


def _first(mask):
    return int(np.argmax(mask)) if mask.any() else None


panels = [
    dict(I=-0.15, title=r'$I=-0.15<I_\ast$', label='A', w0_boost=1.014,
         arrows=[(10, 'theta_gt', -2, 2), (10, 'w_lt', 0.06, 2),
                 (13, 'theta_gt', 2, 1), (18, 'theta_gt', 0, 1),
                 (5, 'theta_gt', -2, 1)],
         static_arrows=[(-2., 0., (1., 0.)), (0., 0., (-1., 0.)), (2., 0., (1., 0.))]),
    dict(I=-0.1069150434, title=r'$I>I_\ast$ but $I \approx I_\ast$', label='B', w0_boost=None,
         arrows=[(10, 'w_lt', 0.05, 2), (10, 'theta_gt', -1.5, 2), (10, 'theta_gt', 3, 2),
                 (5, 'theta_gt', -2, 1), (5, 'w_lt', 0.02, 1),
                 (15, 'theta_gt', 0, 1), (20, 'theta_gt', 0, 1)],
         static_arrows=[(-2., 0., (1., 0.)), (0., 0., (-1., 0.)), (2., 0., (1., 0.))]),
    dict(I=-0.05, title=r'$I_\ast < I= -0.05 < I_c=0$', label='C', w0_boost=None,
         arrows=[(10, 'theta_gt', 0, 2),
                 (5, 'theta_gt', -2, 1), (6, 'theta_gt', 2, 1),
                 (15, 'theta_gt', 0, 1), (20, 'theta_gt', 0, 1)],
         static_arrows='C'),
    dict(I=0.05, title=r'$I=0.05>I_c=0$', label='D', w0_boost=None,
         arrows=[(10, 'theta_gt', 0, 2),
                 (5, 'theta_gt', 0, 1), (15, 'theta_gt', 0, 1), (20, 'theta_gt', 0, 1)],
         static_arrows=None),
]

for p in panels:
    trajs = []
    for ijk in range(1, 26):
        w0 = w_max * ijk / 20 * 2
        if p['w0_boost'] and ijk > 10:
            w0 *= p['w0_boost']
        trajs.append(simulate(p['I'], w0))
    p['trajs'] = trajs

if __name__ == "__main__":

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))

    for ax, p in zip(axes.flat, panels):
        theta_plus, theta_minus = theta_thresholds(p['I'])
        ax.plot(theta_plus, 0, 'ok', markersize=6, markerfacecolor='w')
        ax.plot(theta_minus, 0, 'ok', markersize=6, markerfacecolor='k')

        for ijk, (theta, w) in enumerate(p['trajs'], start=1):
            ax.plot(theta, w, '-k', linewidth=3 if ijk == 10 else 1)

        for ijk, cond, thresh, width in p['arrows']:
            theta, w = p['trajs'][ijk - 1]
            mask = theta > thresh if cond == 'theta_gt' else w < thresh
            ind = _first(mask)
            _arrow(ax, theta, w, ind, width)

        if p['static_arrows'] == 'C':
            static = [(-2., 0., (1., 0.)), (-theta_plus / 3, 0., (-1., 0.)), (2., 0., (1., 0.))]
        else:
            static = p['static_arrows']
        if static:
            for x, y, vec in static:
                draw_arrow(ax, (A, B), (C, D), x, y, np.array(vec), epsilon=0.06, width=1, color='k')

        ax.set_xlim(A, B)
        ax.set_ylim(C, D)
        ax.set_box_aspect(1)
        ax.set_title(p['title'])
        ax.text(-np.pi - 0.9, 0.4, p['label'], fontsize=16, fontweight='bold')

    axes[1, 0].set_xlabel(r'$\theta$')
    axes[1, 0].set_ylabel('$z$')
    axes[1, 1].set_xlabel(r'$\theta$')
    axes[0, 0].set_ylabel('$z$')

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
