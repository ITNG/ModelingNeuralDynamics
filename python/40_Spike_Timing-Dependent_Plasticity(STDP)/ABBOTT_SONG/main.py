import pylab as plt
import numpy as np
from numpy import exp, abs

tau_plus = 10
tau_minus = 10
K_plus = 0.1
K_minus = 2 / 3 * K_plus
z = np.arange(1, 1000) / 1000 * tau_plus * 2


plt.plot(z, K_plus * exp(-z / tau_plus), c='k', lw=3, alpha=0.8)
plt.plot(z, K_plus * exp(-z / tau_plus) * (1 - exp(-abs(z) * 5 / tau_plus)),
         lw=1, c='r')
z = -z
plt.plot(z, -K_minus * exp(z / tau_minus), c='k', lw=3, alpha=0.8)
plt.plot(z, -K_minus * exp(z / tau_minus) * (1 - exp(-abs(z) * 5 / tau_plus)),
         lw=1, c='r')


plt.plot([0, 0], [-K_minus * 1.2, K_plus * 1.2], c='k', ls='--')
plt.plot([-2*tau_plus,2*tau_plus], [0, 0], c='k', ls='--')

plt.xlabel("z [ms]")
plt.title("$F_0$ (black) and $F$ (red)")
plt.margins(x=0, y=0)
plt.show()
