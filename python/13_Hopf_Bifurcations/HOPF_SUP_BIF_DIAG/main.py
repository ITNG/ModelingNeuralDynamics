import numpy as np
import matplotlib.pyplot as plt

I_neg = np.arange(-100, 1) / 100
I_pos = -I_neg

plt.figure(figsize=(6, 6))
plt.plot(I_neg, np.zeros(101), color='k', linewidth=2)
plt.plot(I_pos, np.zeros(101), color='k', linewidth=2, linestyle='dashed')
plt.plot(I_pos, np.sqrt(I_pos), color='k', linewidth=2)
plt.plot(I_pos, -np.sqrt(I_pos), color='k', linewidth=2)

plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.gca().set_box_aspect(1)
plt.xlabel('$I$')
plt.ylabel('oscillation amplitude')
plt.tight_layout()
plt.savefig("fig.png")
# plt.show()
