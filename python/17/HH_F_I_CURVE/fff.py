from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl
from copy import copy
from time import time


i_ext_vec = np.linspace(3, 13, 23)

f_forward = np.loadtxt("forward.txt")

I = f_forward[:, 0]
f = f_forward[:, 1]

pl.plot(I, f, "ko", label="forward")

index = np.where(f == 0)[0]
index = max(index)
I_c = (I[index] + I[index + 1]) / 2.0;
pl.plot([I[index + 1], I[index + 1]], [0, f[index + 1]], '--b', lw=1)



f_backward = np.loadtxt("backward.txt")

I = f_backward[:, 0][::-1]
f = f_backward[:, 1][::-1]

pl.plot(I, f, "ro", fillstyle="none", markersize=8, label="backward")

index = np.where(f == 0)[0]
index = max(index)
I_star = (I[index] + I[index + 1]) / 2.0;

pl.plot([I[index + 1], I[index + 1]], [0, f[index + 1]], '--b', lw=1)
pl.text(I_star - 0.1, -15, r"$I_{\ast}$", fontsize=20, color="b")
pl.text(I_c - 0.1, -15, r"$I_c$", fontsize=20, color="b")

pl.xlim(min(i_ext_vec), max(i_ext_vec))
pl.xlabel("I", fontsize=16)
pl.ylabel("frequency", fontsize=16)
pl.legend(fontsize=14, frameon=False)
pl.tight_layout()
pl.tick_params(labelsize=14)
pl.savefig("fig_17_1.png")