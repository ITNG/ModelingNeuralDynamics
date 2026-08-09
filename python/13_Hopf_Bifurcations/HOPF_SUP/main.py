import numpy as np
import matplotlib.pyplot as plt

r = np.arange(101) / 100 * 1.2

plt.figure(figsize=(7, 7))
for I in [-1, 0, 1]:
    plt.plot(r, I * r - r ** 3, color='k', linewidth=2)

plt.plot([0, 1.2], [0, 0], color='k', linestyle='dashed')
plt.text(0.64, -0.75, '$I=0$', fontsize=16)
plt.text(0.30, -0.95, '$I=-1$', fontsize=16)
plt.text(0.94, -0.55, '$I=1$', fontsize=16)

plt.xlim(0, 1.2)
plt.ylim(-3, 1)
plt.xticks(np.arange(0, 1.21, 0.3))
plt.xlabel('$r$')
plt.ylabel('$f$')
plt.tight_layout()
plt.savefig("fig.png")
# plt.show()
