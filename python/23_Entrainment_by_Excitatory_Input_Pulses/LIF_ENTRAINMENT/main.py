import numpy as np
import matplotlib.pyplot as plt

T = 20.
tau = 30.
n = 10
r = np.exp(-T / tau)
epsilon = 0.4875

t = np.arange(101) / 100 * T
t_final = 1000.

num_pulses = round(t_final / T)

segments = []  # (t_offset, v_post) for each inter-pulse decay curve
reset_lines = []  # t where v hits threshold and resets to 0
dotted_lines = []  # (t, v_pre, v_post) jump at each pulse

v_post = epsilon
for k in range(1, num_pulses + 1):
    segments.append((k * T, v_post))
    v_pre = v_post * r
    v_post = v_pre + epsilon
    v_post = v_post * (v_post < 1)
    if v_post == 0:
        reset_lines.append((k + 1) * T)
    dotted_lines.append(((k + 1) * T, v_pre, v_post))

if __name__ == "__main__":

    plt.figure(figsize=(7, 4))
    plt.plot([0, T], [0, 0], '-k', linewidth=4)
    plt.plot([T, T], [0, epsilon], ':k', linewidth=1)

    for t_offset, v in segments:
        plt.plot(t_offset + t, v * np.exp(-t / tau), '-k', linewidth=2)
    for tt in reset_lines:
        plt.plot([tt, tt], [0, 5], '-k', linewidth=2)
    for tt, v_pre, v_post in dotted_lines:
        plt.plot([tt, tt], [v_pre, v_post], ':k', linewidth=1)

    plt.xlim(0, t_final)
    plt.ylim(0, 6)
    plt.xlabel('$t$ [ms]')
    plt.ylabel('$v$')
    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
