from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

tau_m=10.
I=0.13
w_k=0.05
tilde_w_k=0.08
tau_w=40.

t_final=100.
dt=0.01


def derivative(v,t,wk):
    dv=-v/tau_m+I-wk*np.exp(-t/tau_w)*v
    return dv


t=np.arange(0,t_final+dt,dt)

v=odeint(derivative,0.,t,args=(w_k,))[:,0]
v_tilde=odeint(derivative,0.,t,args=(tilde_w_k,))[:,0]

plt.figure(figsize=(7,4))
plt.plot(t,v,color='k',linewidth=2,label='$v$')
plt.plot(t,v_tilde,color='k',linewidth=2,linestyle='dashed',label=r'$\tilde{v}$')
plt.xlabel('$t$')
plt.title(r'$v$ (solid) and $\tilde{v}$ (dashes)')
plt.xlim(0,60)
plt.ylim(0,1)
plt.tight_layout()
plt.savefig("fig.png")
#plt.show()
