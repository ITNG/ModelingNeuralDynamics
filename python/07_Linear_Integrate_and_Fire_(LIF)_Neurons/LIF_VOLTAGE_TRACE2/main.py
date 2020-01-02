import numpy as np
import pylab as pl


def integrate_rk4(x, dt, f):
    k1 = dt * f(x)
    k2 = dt * f(x + 0.5 * k1)
    k3 = dt * f(x + 0.5 * k2)
    k4 = dt * f(x + k3)

    x = x + (k1 + 2.0 * (k2 + k3) + k4) / 6.0
    return x


def derivative(v):
    dv = -v / tau_m + I
    return dv


t_final = 50.0
tau_m = 2.0
dt = 0.01
I = 1 / (1 - np.exp(-20.0 / tau_m)) / tau_m
print I



if __name__ == "__main__":

    numSteps = int(t_final / dt)
    v = np.zeros(numSteps)
    t = np.arange(0, t_final, dt)

    for i in range(1, numSteps):
        v_new = integrate_rk4(v[i - 1], dt, derivative)

        if v_new <= 1:
            v[i] = v_new
        else:
            v[i] = 0.0


pl.figure(figsize=(7, 3))
pl.xlabel("time [ms]", fontsize=14)
pl.ylabel("v [mV]", fontsize=14)
pl.ylim([0, 2])
pl.xlim(0, max(t))
pl.tight_layout()
pl.plot(t, v, lw=2, c="k")
pl.savefig("fig_7_5.png")
pl.show()