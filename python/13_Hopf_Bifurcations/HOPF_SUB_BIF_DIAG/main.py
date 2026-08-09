import numpy as np
import matplotlib.pyplot as plt

I = np.arange(-100, 1) / 100

plt.figure(figsize=(6, 6))
plt.plot(I, np.zeros(101), color='k', linewidth=2)
plt.plot(-I, np.zeros(101), color='k', linewidth=2, linestyle='dashed')
plt.plot(I, np.sqrt(-I), color='k', linewidth=2, linestyle='dashed')
plt.plot(I, -np.sqrt(-I), color='k', linewidth=2, linestyle='dashed')

plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.gca().set_box_aspect(1)
plt.xlabel('$I$')
plt.ylabel('oscillation amplitude')
plt.tight_layout()
plt.savefig("fig.png")
# plt.show()
