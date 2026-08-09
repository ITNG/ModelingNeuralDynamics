import numpy as np
import matplotlib.pyplot as plt

v=np.arange(-1000,501)/10.
r=(120-v)/(1+np.exp(-(v+15)/5))*4/25

plt.figure(figsize=(7,4))
plt.plot(v,r,color='k',linewidth=2)
plt.xlabel('$v$ [mV]')
plt.ylabel(r'$c_\infty$')
plt.tight_layout()
plt.savefig("fig.png")
#plt.show()
