import numpy as np
import matplotlib.pyplot as plt

r = np.arange(101) / 100 * 1.2

plt.figure(figsize=(7, 7))
for I in [-0.2, 0, 0.2, -0.4]:
    plt.plot(r, I * r + r ** 3 - r ** 5, color='k', linewidth=2)

plt.plot([0, 1.2], [0, 0], color='k', linestyle='dashed')
plt.text(0.7, 0.23, '$I=0$', fontsize=16)
plt.text(0.6, 0.07, '$I=-0.2$', fontsize=16)
plt.text(0.8, 0.39, '$I=0.2$', fontsize=16)
plt.text(0.5, -0.06, '$I=-0.4$', fontsize=16)

plt.xlim(0, 1.2)
plt.ylim(-0.5, 0.5)
plt.xticks(np.arange(0, 1.21, 0.3))
plt.xlabel('$r$')
plt.ylabel('$f$')
plt.tight_layout()
plt.savefig("fig.png")
# plt.show()
