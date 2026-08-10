import matplotlib.pyplot as plt

# schematic recreation of a hand-drawn matlab figure (figure.graffle):
# a train of narrow square current pulses of height i_h above a baseline
# i_l, width epsilon, and period T.

T = 3.
epsilon = 0.4
i_l = 1.
i_h = 3.

if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(8, 5))

    t = [0.]
    i = [i_l]
    for k in range(3):
        t0 = k * T
        t += [t0, t0, t0 + epsilon, t0 + epsilon]
        i += [i_l, i_h, i_h, i_l]
    t.append(3 * T - epsilon)
    i.append(i_l)

    ax.plot(t, i, '-k', linewidth=3)
    ax.annotate('', xy=(3 * T, 0), xytext=(-0.2, 0),
                arrowprops=dict(arrowstyle='->', linewidth=2))
    ax.annotate('', xy=(0, i_h + 1), xytext=(0, -0.1),
                arrowprops=dict(arrowstyle='->', linewidth=2))

    ax.set_xlim(-0.3, 3 * T)
    ax.set_ylim(-0.7, i_h + 1.2)
    ax.axis('off')

    xticks = [-0.15, epsilon + 0.15, T - 0.15, T + epsilon + 0.25, 2 * T - 0.15, 2 * T + epsilon + 0.25]
    xticklabels = ['0', r'$\epsilon$', '$T$', r'$T\!+\!\epsilon$', '$2T$', r'$2T\!+\!\epsilon$']
    for xt, label in zip(xticks, xticklabels):
        ax.text(xt, -0.6, label, ha='center', fontsize=14)
    for yt, label in zip([i_l, i_h], ['$I_L$', '$I_H$']):
        ax.text(-0.35, yt, label, ha='right', va='center', fontsize=16)

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
