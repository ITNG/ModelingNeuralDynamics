import numpy as np
import matplotlib.pyplot as plt

tau_m=10.
I=0.12
tau_w=100.
delta=0.01
z_max=0.05
period=tau_m*np.log(tau_m*I/(tau_m*I-1))
frequency=1000/period
print(frequency)

dt=0.01
dt05=dt/2

N=100

v=np.zeros(N+1)
w=np.arange(N+1)/N*z_max
phi=np.zeros(N+1)
t=0.

done=np.zeros(N+1,dtype=bool)

while not done.all():
    v_old=v.copy()
    w_old=w.copy()
    v_inc=-v/tau_m+I-w*v
    w_inc=-w/tau_w
    v_tmp=v+dt05*v_inc
    w_tmp=w+dt05*w_inc
    v_inc=-v_tmp/tau_m+I-w_tmp*v_tmp
    w_inc=-w_tmp/tau_w
    v=v+dt*v_inc
    w=w+dt*w_inc
    t=t+dt

    ind=(v>1)&(~done)
    done[ind]=True
    phi[ind]=(v[ind]-1)*w_old[ind]+(1-v_old[ind])*w[ind]
    phi[ind]=phi[ind]/(v[ind]-v_old[ind])+delta

z=np.arange(N+1)/N*z_max

plt.figure(figsize=(6,6))
plt.plot(z,phi,color='k',linewidth=2)
plt.plot([0,z_max],[0,z_max],color='k',linestyle='dashed')

z_star=0.02791
plt.plot([z_star,z_star],[0,z_star],color='k',linestyle='dotted',linewidth=1)
plt.plot(z_star,0,'.r',markersize=20)
plt.text(0.0265,-0.003,r'$z_\ast$',fontsize=16,color='r')

plt.xlim(0,z_max)
plt.ylim(0,z_max)
plt.gca().set_aspect('equal')
plt.xticks([0,0.02,0.04])
plt.yticks([0,0.02,0.04])
plt.xlabel('$z$')
plt.ylabel(r'$\phi(z)$')
plt.tight_layout()
plt.savefig("fig.png")
#plt.show()
