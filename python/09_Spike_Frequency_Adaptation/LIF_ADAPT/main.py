import numpy as np
import matplotlib.pyplot as plt

tau_m=10.
I=0.13
tau_a=40.
delta=0.05

t_final=300.
dt=0.01


def derivative(state):
    v,w=state
    dv=-v/tau_m+I-w*v
    dw=-w/tau_a
    return np.array([dv,dw])


def integrate_rk4(x,dt,f):
    k1=dt*f(x)
    k2=dt*f(x+0.5*k1)
    k3=dt*f(x+0.5*k2)
    k4=dt*f(x+k3)
    return x+(k1+2.*(k2+k3)+k4)/6.


num_steps=int(t_final/dt)
t=np.arange(0,t_final,dt)
v=np.zeros(num_steps)
w=np.zeros(num_steps)

for i in range(1,num_steps):
    v_new,w_new=integrate_rk4(np.array([v[i-1],w[i-1]]),dt,derivative)
    if v_new<=1:
        v[i]=v_new
        w[i]=w_new
    else:
        v[i]=0.
        w[i]=w_new+delta

fig,ax=plt.subplots(2,1,figsize=(8,6),sharex=True)
ax[0].plot(t,v,color='k',linewidth=2)
ax[0].set_ylabel('$v$')
ax[0].set_ylim(0,4)
ax[0].set_yticks([0,1])

ax[1].plot(t,w,color='k',linewidth=2)
ax[1].set_xlabel('$t$')
ax[1].set_ylabel('$w$')
ax[1].set_xlim(0,t_final)
ax[1].set_ylim(0,0.1)

plt.tight_layout()
plt.savefig("fig.png")
#plt.show()
