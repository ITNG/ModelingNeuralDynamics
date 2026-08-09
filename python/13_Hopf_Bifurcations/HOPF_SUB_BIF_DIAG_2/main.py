import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(6, 6))

# I in [-1, -0.25]: only the r=0 branch (stable)
I = -1 + np.arange(101) / 100 * 0.75
plt.plot(I, np.zeros(101), color='k', linewidth=2)

# I in [-0.25, 0]: r=0 stable, plus an unstable (r0) and a stable (R0)
# limit-cycle branch
I = -0.25 + np.arange(101) / 100 * 0.25
r0 = np.sqrt(1 / 2 - np.sqrt(1 / 4 + I))
R0 = np.sqrt(1 / 2 + np.sqrt(1 / 4 + I))
plt.plot(I, np.zeros(101), color='k', linewidth=2)
plt.plot(I, r0, color='k', linewidth=2, linestyle='dashed')
plt.plot(I, -r0, color='k', linewidth=2, linestyle='dashed')
plt.plot(I, R0, color='k', linewidth=2)
plt.plot(I, -R0, color='k', linewidth=2)

# I in [0, 1]: r=0 now unstable, only the stable R0 branch remains
I = np.arange(101) / 100
R0 = np.sqrt(1 / 2 + np.sqrt(1 / 4 + I))
plt.plot(I, np.zeros(101), color='k', linewidth=2, linestyle='dashed')
plt.plot(I, R0, color='k', linewidth=2)
plt.plot(I, -R0, color='k', linewidth=2)

plt.xlim(-1, 1)
plt.ylim(-2, 2)
plt.gca().set_box_aspect(1)
plt.xlabel('$I$')
plt.ylabel('oscillation amplitude')
plt.tight_layout()
plt.savefig("fig.png")
# plt.show()
