from scipy.integrate import odeint
from math import floor
from numpy import exp
from copy import copy
import numpy as np
import pylab as pl

# ------------------------------------------------------------------#


def alpha_h(v):
    return (0.128 * exp(-(v + 50.0) / 18.0))


def alpha_m(v):
    return (0.32 * (v + 54.0) / (1.0 - exp(-(v + 54.0) / 4.0)))


def alpha_n(v):
    return (0.032 * (v + 52) / (1.0 - exp(-(v + 52.0) / 5.0)))


def beta_h(v):
    return (4.0 / (1.0 + exp(-(v + 27.0) / 5.0)))


def beta_m(v):
    return (0.28 * (v + 27.0) / (exp((v + 27.0) / 5.0) - 1.0))


def beta_n(v):
    return (0.5 * exp(-(v + 57.0) / 40.0))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))
# ------------------------------------------------------------------#


def derivative(x0, t=0, n0=1):
    '''
    define Traub Model
    '''

    if n0 == 1:
        v, m, n, h, = x0
    else:
        v, m = x0[:n0], x0[n0: (2 * n0)]
        n, h = x0[(2 * n0): (3 * n0)], x0[(3 * n0):]

    dv = i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    if n0 == 1:
        return np.array([dv, dm, dn, dh])
    else:
        return np.hstack((dv, dm, dn, dh))
# ------------------------------------------------------------------#


def rungeKuttaIntegrator(x, dt, f):

    k1 = dt * f(x)
    k2 = dt * f(x + 0.5 * k1)
    k3 = dt * f(x + 0.5 * k2)
    k4 = dt * f(x + k3)

    x = x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    return x
# ------------------------------------------------------------------#


def splayState(i_ext, phiVec, x0, N, f):
    """
    input: i_ext=external drive (a single number)
        phi_vec = vector of phases at which
                    neurons are to be initialized
        The length, num, of i_ext is the total number
        of neurons.
    output: a N-by-3 array called rtm_init. The columns
            contain values of v, h, and n. If i_ext is below
            the firing threshold, then the points (v,h,n)
            (the rows of the matrix rtm_init) are all equal
            to each other, and equal to the stable equilibrium
            point.
            If i_ext is above the firing threshold, then the
            (v,h,n) lie on the limit cycle, with the phases
            given by phi_vec.
            The function also computes the period, T.
    """

    t = 0.0
    dt = 0.01
    t_final = 5000.0
    numSpikes = 0
    iteration = 0
    tSpikes = []
    numSteps = int(t_final / dt)

    v, m = np.zeros(numSteps), np.zeros(numSteps)
    n, h = np.zeros(numSteps), np.zeros(numSteps)

    # ofile = open("v.txt", "w")
    v[0], m[0], n[0], h[0] = x0
    i = 1
    while (numSpikes < 5) and (t < t_final):
        vpre, mpre, npre, hpre = v, m, n, h
        v[i], m[i], n[i], h[i] = rungeKuttaIntegrator(x0, dt, f)

        # ofile.write("%18.4f %18.4f \n" % (i*dt, v))

        x0 = np.array([v[i], m[i], n[i], h[i]])

        if (v[i-1] < spikeThreshold) & (v[i] >= spikeThreshold):
            numSpikes += 1
            tmp = ((i - 1) * dt * (v[i-1] - spikeThreshold) +
                   i * dt * (spikeThreshold - v[i])) / (v[i-1] - v[i])
            tSpikes.append(tmp)

        i += 1
        t = (i + 1) * dt

    # ofile.close()
    # r = np.loadtxt("v.txt")
    # pl.plot(r[:, 0], r[:, 1], lw=1, color="k")
    # pl.show()

    initialConition = np.zeros((N, 3))

    if numSpikes < 5:
        initialConition[:, 0] = v[i-1]
        initialConition[:, 1] = h[i-1]
        initialConition[:, 2] = n[i-1]
        T = np.inf

    elif (numSpikes == 5):
        T = tSpikes[-1] - tSpikes[-2]
        for i in range(N):
            phi0 = phiVec[i]
            t0 = phi0 * T + tSpikes[-2]
            k = int(floor(t0 / dt))
            initialConition[i, 0] = (
                v[k + 1] * (t0 - (k - 1) * dt) + v[k] * (k * dt - t0)) / dt
            initialConition[i, 1] = (
                n[k + 1] * (t0 - (k - 1) * dt) + n[k] * (k * dt - t0)) / dt
            initialConition[i, 2] = (
                h[k + 1] * (t0 - (k - 1) * dt) + h[k] * (k * dt - t0)) / dt

    print ("Period is %10.3f [ms]" % T)

    return initialConition
# ------------------------------------------------------------------#


N = 30
c = 1
g_k = 80
g_na = 100
g_l = 0.1
v_k = -100
v_na = 50
v_l = -67
i_ext = 0.3
t_final = 200.0
dt = 0.01
spikeThreshold = -20.0

v = -70.0
m = m_inf(v)
n = n_inf(v)
h = h_inf(v)
x0 = np.array([v, m, n, h])

if __name__ == "__main__":

    phiVec = np.arange(N - 1.0, -1.0, -1.0) / float(N) + 1.0 / (2.0 * N)
    initialVec = splayState(i_ext,
                            phiVec,
                            x0,
                            N,
                            derivative)

    v = initialVec[:, 0]
    m = m_inf(v)
    n = initialVec[:, 1]
    h = initialVec[:, 2]

    x0 = np.hstack((v, m, n, h))
    t = np.arange(0, t_final, 0.01)
    sol = odeint(derivative, x0, t, args=(N,))
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
        ax.plot(tSpikes[i], [i + 1] * len(tSpikes[i]), "k.")

    ax.set_xlim(0, t_final)
    ax.set_xlabel("time [ms]", fontsize=14)
    ax.set_ylabel("neuron #", fontsize=14)
    pl.tight_layout()
    pl.savefig("fig-24-1.png")
    pl.close()
