import numpy as np
import matplotlib.pyplot as plt

tau_m=2.
I=0.15

t_final=150
dt=0.01
dt05=dt/2
num_steps=int(t_final/dt)
t=np.linspace(0.0,t_final,num_steps)
v=np.zeros_like(t)
plt.figure(figsize=(10,5))
for i in range(num_steps-1): 
    v_inc=-v[i]/tau_m*(1.0-v[i])+I;
    v_tmp=v[i]+dt05*v_inc;
    v_inc=-v_tmp/tau_m*(1-v_tmp)+I;
    v_new=v[i]+dt*v_inc; 
    if v_new<=1.0:
        v[i+1]=v_new
        plt.plot([t[i],t[i+1]],[v[i],v[i+1]],color='black',linestyle='solid')
    else:
        v[i+1]=0.0
        plt.plot([t[i],t[i+1]],[v[i],v[i+1]],color='black',linestyle='dashed')
#plt.plot(t,v)
plt.xlabel('$t$')
plt.ylabel('$v$')
plt.ylim((0.0,2.0))
plt.xlim((0.0,150.0))
plt.savefig("fig.png")
#plt.show()
exit()
