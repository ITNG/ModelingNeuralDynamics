import numpy as np
import matplotlib.pyplot as plt

tau_m=1./2.
I=0.505
t_final=150

dt=0.001
dt05=dt/2.0

m_steps=int(t_final/dt)
t=np.linspace(0,t_final,m_steps+1)
theta=np.zeros(m_steps+1)

for k in range(m_steps):
    theta_inc=-np.cos(theta[k])/tau_m+2.*I*(1.+np.cos(theta[k]))
    theta_tmp=theta[k]+dt05*theta_inc
    theta_inc=-np.cos(theta_tmp)/tau_m+2.*I*(1.+np.cos(theta_tmp))
    theta[k+1]=theta[k]+dt*theta_inc


plt.figure(figsize=(10,5))
plt.plot(t,1.0-np.cos(theta))
plt.xlabel(r'$t$')
plt.ylabel(r'$1-\cos( \theta )$')
plt.ylim((0.0,2.0))
plt.xlim((0.0,150.0))
plt.xticks(np.arange(0, 151., step=50))
plt.yticks(np.arange(0, 2.1, step=.50))
plt.grid()
plt.savefig("fig.png")
#plt.show()
