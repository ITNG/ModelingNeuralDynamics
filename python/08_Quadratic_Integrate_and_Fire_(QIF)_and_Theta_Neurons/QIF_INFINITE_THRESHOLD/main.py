import numpy as np
import matplotlib.pyplot as plt

tau_m=2.
I=0.15

w=np.sqrt(tau_m*I-1./4.)
T=2.*tau_m/w*np.arctan(1./(2.*w))
t_ast=tau_m/w*(np.pi/2.-np.arctan(1./(2.*w)))

t_period=np.arange(101)/100.*T
v_0_to_1=0.5+w*np.tan(w/tau_m*t_period-np.arctan(1./(2.*w)))

t_blowup=np.arange(100)/100.*t_ast
v_1_to_inf=0.5+w*np.tan(w/tau_m*t_blowup+np.arctan(1./(2.*w)))

v_minus_inf_to_0=1.-v_1_to_inf[::-1]

plt.figure(figsize=(10,5))
plt.plot(t_period,v_0_to_1,color='black',linewidth=3)
plt.plot(T+t_blowup,v_1_to_inf,color='black',linewidth=1)
plt.plot([T+t_ast,T+t_ast],[-1,2],color='black',linestyle='dashed',linewidth=1)

for ijk in range(1,6):
    plt.plot(ijk*(T+2*t_ast)-t_ast+t_blowup,v_minus_inf_to_0,color='black',linewidth=1)
    plt.plot(ijk*(T+2*t_ast)+t_period,v_0_to_1,color='black',linewidth=3)
    plt.plot(ijk*(T+2*t_ast)+T+t_blowup,v_1_to_inf,color='black',linewidth=1)
    plt.plot([ijk*(T+2*t_ast)+T+t_ast,ijk*(T+2*t_ast)+T+t_ast],[-1,2],color='black',linestyle='dashed',linewidth=1)

plt.xlabel('$t$')
plt.ylabel('$v$')
plt.ylim((-1,2))
plt.xlim((0,150))
plt.savefig("fig.png")
#plt.show()
