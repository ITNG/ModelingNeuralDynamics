import numpy as np
import matplotlib.pyplot as plt

from mnd.core import draw_arrow

theta=np.linspace(0,2*np.pi,101)
x=np.cos(theta)
y=np.sin(theta)

fig,axes=plt.subplots(1,3,figsize=(12,4.5))

# panel 1: I < 1/(4 tau_m) -- one stable (black) and one unstable (white) fixed point
ax=axes[0]
ax.plot(x,y,color='k',linewidth=2)
eps=0.15
theta0=-0.4*np.pi
x0,y0=np.cos(theta0),np.sin(theta0)
ax.fill(x0+eps*x,y0+eps*y,'k')
y0=-y0
ax.fill(x0+eps*x,y0+eps*y,'w',linewidth=1,edgecolor='k')

for theta0 in [0.4,-0.4]:
    x0,y0=np.cos(theta0),np.sin(theta0)
    v=np.array([np.sin(theta0+0.1),-np.cos(theta0+0.1)])
    draw_arrow(ax,[-1.5,1.5],[-1.5,1.5],x0,y0,v)

for theta0 in [-2/3*np.pi,2/3*np.pi,np.pi]:
    x0,y0=np.cos(theta0),np.sin(theta0)
    v=-np.array([np.sin(theta0-0.1),-np.cos(theta0-0.1)])
    draw_arrow(ax,[-1.5,1.5],[-1.5,1.5],x0,y0,v)

ax.text(-1.0,1.6,r'$I < 1/(4\tau_m)$',fontsize=14)
ax.axis([-1.5,1.5,-1.5,1.5])
ax.set_aspect('equal')
ax.axis('off')

# panel 2: I = 1/(4 tau_m) -- a single half-stable fixed point
ax=axes[1]
ax.plot(x,y,color='k',linewidth=2)
theta_half=np.linspace(0,np.pi,101)
ax.fill(np.concatenate(([1-eps,1+eps],1+eps*np.cos(theta_half))),
        np.concatenate(([0,0],-eps*np.sin(theta_half))),'k')
ax.fill(np.concatenate(([1-eps,1+eps],1+eps*np.cos(theta_half))),
        np.concatenate(([0,0],eps*np.sin(theta_half))),'w',linewidth=1,edgecolor='k')

for theta0 in [-2/3*np.pi,1/3*np.pi,-1/3*np.pi,2/3*np.pi,np.pi]:
    x0,y0=np.cos(theta0),np.sin(theta0)
    v=-np.array([np.sin(theta0-0.1),-np.cos(theta0-0.1)])
    draw_arrow(ax,[-1.5,1.5],[-1.5,1.5],x0,y0,v)

ax.text(-1.0,1.6,r'$I = 1/(4\tau_m)$',fontsize=14)
ax.axis([-1.5,1.5,-1.5,1.5])
ax.set_aspect('equal')
ax.axis('off')

# panel 3: I > 1/(4 tau_m) -- no fixed points, the neuron fires
ax=axes[2]
ax.plot(x,y,color='k',linewidth=2)

for theta0 in [-2/3*np.pi,1/3*np.pi,-1/3*np.pi,2/3*np.pi,np.pi,0]:
    x0,y0=np.cos(theta0),np.sin(theta0)
    v=-np.array([np.sin(theta0-0.1),-np.cos(theta0-0.1)])
    draw_arrow(ax,[-1.5,1.5],[-1.5,1.5],x0,y0,v)

ax.text(-1.1,1.6,r'$I > 1/(4\tau_m)$',fontsize=14)
ax.axis([-1.5,1.5,-1.5,1.5])
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.savefig("fig.png")
#plt.show()
