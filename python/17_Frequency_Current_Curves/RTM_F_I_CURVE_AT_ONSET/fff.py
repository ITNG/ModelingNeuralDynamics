from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl
from copy import copy
from time import time


i_ext_high = 0.1194
i_ext_low = 0.1193
i_ext_vec = np.linspace(i_ext_low, i_ext_high, 11)


f_forward = np.loadtxt("forward.txt")
I = f_forward[:, 0]
f = f_forward[:, 1]

pl.plot(I, f, "ko", label="forward")


# plot the solid red line
ind = np.where(f > 0)[0]
I0 = I[ind]
f0 = f[ind]
alpha_vec = np.linspace(0, 1, 101)

I_c_low = i_ext_vec[min(ind)-1]
I_c_high = i_ext_vec[min(ind)]

C_vec = np.zeros(len(alpha_vec))
err_vec = np.zeros(len(alpha_vec))

for ijk in range(len(alpha_vec)):

    alpha = alpha_vec[ijk]
    I_c = I_c_low * alpha + I_c_high * (1 - alpha)
    y = f0 / np.sqrt(I0 - I_c)
    print (y)
    exit(0)
    # measure how constant f/sqrt(I-I_c) is, with
    # the I_c defined above.
    err_vec[ijk] = (max(y) - min(y)) / np.mean(y)
    C_vec[ijk] = np.mean(y)

# pick the I_c that makes f/sqrt(I-I_c) as constant
ind = np.argmin(err_vec)
alpha = alpha_vec[ind]

I_c = I_c_low * alpha + I_c_high * (1 - alpha)
C = C_vec[ind]
I = I_c + np.linspace(0, 1.0, 1001) * (i_ext_high - I_c)
print (C * np.sqrt(I - I_c))
pl.plot(I, C * np.sqrt(I - I_c), 'r', lw=2)


f_backward = np.loadtxt("backward.txt")
I = f_backward[:, 0][::-1]
f = f_backward[:, 1][::-1]

pl.plot(I, f, "ro", fillstyle="none", markersize=8, label="backward")

# pl.xlim(min(i_ext_vec), max(i_ext_vec))
pl.xlabel("I", fontsize=16)
pl.ylabel("frequency", fontsize=16)
pl.legend(fontsize=14, frameon=False)
pl.tight_layout()
pl.tick_params(labelsize=14)
pl.savefig("fig_17_5.png")
