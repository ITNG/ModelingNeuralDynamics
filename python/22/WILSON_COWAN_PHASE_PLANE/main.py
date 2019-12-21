import numpy as np
import pylab as pl
from copy import copy
from scipy.integrate import odeint


def g(x):
    return 100 * x * x / (400.0 + x * x) * (x > 0)


def f(x):
    return 100.0 * x * x / (900 + x * x) * (x > 0)


def dE(E, I):
    return (f(w_EE * E - w_IE * I + I_E) - E) / tau_E


def dI(E, I):
    return (g(w_EI * E - w_II * I + I_I) - I) / tau_I


def derivative(x, t):

    E, I = x
    return [dE(E, I), dI(E, I)]


def plot_streamplot(ax,
                    derivative,
                    x=[0, 100],
                    y=[0, 100],
                    dx=0.1,
                    dy=0.1):

    phi1, phi2 = np.meshgrid(np.arange(x[0], x[1], dx),
                             np.arange(y[0], y[1], dy))
    dphi1_dt, dphi2_dt = derivative([phi1, phi2], 0)
    ax.streamplot(phi1, phi2,
                  dphi1_dt, dphi2_dt,
                  color='k',
                  linewidth=0.5,
                  cmap=pl.cm.autumn)
    # ax.quiver(phi1, phi2, dphi1_dt, dphi2_dt)


def inv_f(x):
    return np.sqrt(900 * x / (100.0 - x)) * (x < 100.0)


def inv_g(x):
    return np.sqrt(400 * x / (100.0 - x)) * (x < 100.0)

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

    E = [0, 100]
    I = [0, 100]

    fig, ax = pl.subplots(figsize=(5, 5))
    plot_streamplot(ax, derivative, E, I, 0.1, 0.1)

    
    def plot_nullcline(f, ax, label=None, color="r"):
        x = 100.0
        I_left = np.arange(0, x, 0.1)
        I_right = np.arange(0.1, x + 0.1, 0.1)
        nx = len(I_left)
        i_red_index = 0
        E_red = []
        I_red = []
        # bisection method to find zero
        for i in range(nx):
            E = i / float(nx) * 100.0
            R_left = f(E, I_left)
            R_right = f(E, I_right)
            ind = np.where((R_left * R_right) < 0)[0]
            if len(ind) > 0:
                for j in range(len(ind)):
                    I_l = I_left[ind[j]]
                    I_r = I_right[ind[j]]
                    while (I_r-I_l) > 1e-8:
                        I_c = 0.5 * (I_r + I_l)
                        R_l = f(E, I_l)
                        R_c = f(E, I_c)

                        if (R_l * R_c) < 0:
                            I_r = I_c
                        else:
                            I_l = I_c

                        I_c = 0.5 * (I_r + I_l)
                        i_red_index = i_red_index + 1.0
                        E_red.append(E)
                        I_red.append(I_c)

        ax.plot(E_red, I_red, lw=2, c=color, ls="--", label=label)
        if label:
            ax.legend()
    
    plot_nullcline(dE, ax, label="dE/dt=0", color="r")
    plot_nullcline(dI, ax, label="dI/dt=0", color="b")

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    ax.set_xlabel("E", fontsize=14)
    ax.set_ylabel("I", fontsize=14)
    ax.tick_params(labelsize=14)
    pl.tight_layout()
    pl.savefig("fig_22_3.png")
    # pl.show()
