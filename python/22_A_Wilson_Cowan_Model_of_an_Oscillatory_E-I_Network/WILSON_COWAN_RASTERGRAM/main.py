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

def extract_spikes(x, t, num_nodes):

        spikes = []
        num_steps = len(x)

        for m in range(num_nodes):
            spk = []
            for i in range(num_steps - 1):
                r = np.random.rand()
                if r < 5e-4 * (x[i] + x[i + 1]) * dt:
                    spk.append(t[i])
            spikes.append(spk)

        return spikes



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

    num_i = 20
    num_e = 80
    
    # extract spikes from activities
    spikes_e = extract_spikes(E, t, num_e)
    spikes_i = extract_spikes(I, t, num_i)
    
    fig , ax = pl.subplots(1, figsize=(7, 4))

    # plot spikes
    for ii in range(num_e):
        ax.plot(spikes_e[ii], [ii+num_i] *
                len(spikes_e[ii]), '.',
                c='r',
                markersize=2)

    for ii in range(num_i):
        ax.plot(spikes_i[ii], [ii] *
                len(spikes_i[ii]), '.',
                c='b',
                markersize=2)

    ax.set_xlim(min(t), max(t))
    ax.set_ylim(0, num_e+num_i)
    ax.set_xlabel("time [ms]", fontsize=14)
    ax.set_ylabel("v [mV]", fontsize=14)
    # pl.yticks([0, 50, 100])
    pl.tick_params(labelsize=14)
    pl.tight_layout()
    pl.savefig("fig_22_4.png")
    # pl.show()
