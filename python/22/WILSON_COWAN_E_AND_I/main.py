import numpy as np
import pylab as pl
from scipy.integrate import odeint


def g(x):
    if x >= 0:
        return 100 * x * x / (400.0 + x * x)
    else:
        return 0.0

def f(x):
    if x >= 0:
        return 100.0 * x * x / (900 + x * x)
    else:
        return 0.0

def derivative(x, t):

    E, I = x

    dE = (f(w_EE * E - w_IE * I + I_E) - E) / tau_E
    dI = (g(w_EI * E - w_II * I + I_I) - I) / tau_I

    return [dE, dI]


I_E = 20.0
I_I = 0.0
w_EE = 1.5
w_IE = 1.0
w_EI = 1.0
w_II = 0.0
tau_E = 5.0
tau_I = 10.0

t_final = 300.0
dt = 0.01

E0 = 50.0
I0 = 10.0


if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, [E0, I0], t)
    E = sol[:, 0]
    I = sol[:, 1]

    pl.figure(figsize=(7, 3))
    pl.plot(t, E, lw=2, c="r", label="E")
    pl.plot(t, I, lw=2, c="b", label="I")
    pl.xlim(min(t), max(t))
    pl.ylim(0, 100)
    pl.xlabel("time [ms]", fontsize=14)
    pl.ylabel("v [mV]", fontsize=14)
    pl.yticks([0, 50, 100])
    pl.tick_params(labelsize=14)
    pl.tight_layout()
    pl.legend(fontsize=14)
    pl.savefig("fig_22_2.png")
    # pl.show()
