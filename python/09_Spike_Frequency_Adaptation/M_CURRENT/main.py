import numpy as np
import matplotlib.pyplot as plt

v=np.arange(-100,51)
w_inf=1./(1+np.exp(-(v+35)/10))
tau_w=400./(3.3*np.exp((v+35)/20)+np.exp(-(v+35)/20))

fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].plot(v,w_inf,color='k',linewidth=2)
ax[0].set_xlabel('$v$ [mV]')
ax[0].set_ylabel(r'$w_\infty$')

ax[1].plot(v,tau_w,color='k',linewidth=2)
ax[1].set_xlabel('$v$ [mV]')
ax[1].set_ylabel(r'$\tau_w$')

plt.tight_layout()
plt.savefig("fig.png")
#plt.show()
