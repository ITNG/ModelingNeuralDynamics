from scipy.integrate import odeint
from math import floor
from numpy import exp
from copy import copy
import numpy as np
import pylab as pl
import lib


# ------------------------------------------------------------------#
N = 30
c = 1.0
g_k = 80.0
g_na = 100.0
g_l = 0.1
v_k = -100.0
v_na = 50.0
v_l = -67.0
i_ext = 0.3
t_final = 200.0
dt = 0.01
spikeThreshold = -20.0

tau_d = 2.0
tau_r = 0.5
tau_peak = 0.5
g_syn = 0.0075


if __name__ == "__main__":
    v = -70.0
    m = lib.m_inf(v)
    n = lib.n_inf(v)
    h = lib.h_inf(v)

    x0 = np.array([v, m, n, h])
    phiVec = np.arange(N - 1.0, -1.0, -1.0) / float(N) + 1.0 / (2.0 * N)
    initialVec = lib.splayState(i_ext,
                                phiVec,
                                x0,
                                N,
                                lib.derivativeSingle)

    v = initialVec[:, 0]
    m = lib.m_inf(v)
    n = initialVec[:, 1]
    h = initialVec[:, 2]
    q = np.zeros(N)
    s = np.zeros(N)
    x0 = np.hstack((v, m, n, h, q, s))

    t = np.arange(0, t_final, 0.01)
    sol = odeint(lib.derivativePopulation, x0, t, args=(N,))
    v = sol[:, :N]

    tSpikes = []
    nSteps = len(t)
    for i in range(N):
        tspk = []
        for j in range(1, nSteps - 1):
            if (v[j-1][i] < spikeThreshold) & (v[j][i] >= spikeThreshold):

                ts = ((j - 1) * dt * (v[j-1][i] - spikeThreshold) +
                      j * dt * (spikeThreshold - v[j][i])) / (v[j-1][i] - v[j][i])
                tspk.append(ts)
        tSpikes.append(tspk)

    fig, ax = pl.subplots(1, figsize=(10, 4))
    for i in range(N):
        ax.plot(tSpikes[i], [i + 1] * len(tSpikes[i]), "r.")

    ax.set_xlim(0, t_final)
    ax.set_xlabel("time [ms]", fontsize=14)
    ax.set_ylabel("neuron #", fontsize=14)
    pl.tight_layout()
    pl.savefig("fig-24-3.png")
    pl.close()
