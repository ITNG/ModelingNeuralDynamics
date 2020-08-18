import numpy as np
from numpy import exp
import pylab as pl


def tau_peak_function(tau_d, tau_r, tau_d_q, dt):

    # dt = 0.01
    dt05 = 0.5 * dt

    s = 0
    t = 0
    s_inc = exp(-t / tau_d_q) * (1.0 - s) / tau_r - s * tau_d
    while s_inc > 0:
        t_old = t
        s_inc_old = s_inc
        s_tmp = s + dt05 * s_inc
        s_inc_tmp = exp(-(t + dt05) / tau_d_q) * \
            (1.0 - s_tmp) / tau_r - s_tmp / tau_d
        s = s + dt * s_inc_tmp
        t = t + dt
        s_inc = exp(-t / tau_d_q) * (1.0 - s) / tau_r - s / tau_d

    return (t_old * (-s_inc) + t * s_inc_old) / (s_inc_old - s_inc)


def tau_d_q_function(tau_d, tau_r, tau_hat, dt):

    # set an interval for tau_d_q
    tau_d_q_left = 1.0
    while tau_peak_function(tau_d, tau_r, tau_d_q_left, dt) > tau_hat:
        tau_d_q_left *= 0.5

    tau_d_q_right = tau_r
    while tau_peak_function(tau_d, tau_r, tau_d_q_right, dt) < tau_hat:
        tau_d_q_right *= 2.0

    # bisection method
    while tau_d_q_right - tau_d_q_left > 1e-12:
        tau_d_q_mid = 0.5 * (tau_d_q_left + tau_d_q_right)
        if (tau_peak_function(tau_d, tau_r, tau_d_q_mid, dt) <= tau_hat):
            tau_d_q_left = tau_d_q_mid
        else:
            tau_d_q_right = tau_d_q_mid

    return 0.5 * (tau_d_q_left + tau_d_q_right)


def splayState(i_ext, phiVec, x0, f, N=1):
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
