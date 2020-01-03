from scipy.integrate import odeint
from math import floor
from numpy import exp
from copy import copy
import numpy as np
import pylab as pl
from main import *


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


def derivativeSingle(x0, t=0):
    '''
    define Traub Model for a single neuron
    '''

    v, m, n, h, = x0

    dv = i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return np.array([dv, dm, dn, dh])


def derivativePopulation(x0, t=0, n0=2):
    '''
    define Traub Model for a population of neurons
    '''

    v, m = x0[:n0], x0[n0: (2 * n0)]
    n, h = x0[(2 * n0) : (3 * n0)], x0[(3 * n0) : (4 * n0)]
    q, s = x0[(4 * n0): (5 * n0)], x0[(5 * n0): (6 * n0)]

    I_Na = g_na * h * m ** 3 * (v - v_na)
    I_K = g_k * n ** 4 * (v - v_k)
    I_L = g_l * (v - v_l)
    I_syn = g_syn * (s - np.sum(s)) * v
    
    dv = i_ext - I_Na - I_K - I_L + I_syn
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h
    dq = 5.0 * (1.0 + np.tanh(0.1 * v)) * (1.0 - q) - q / tau_d_q
    ds = q * (1 - s) / tau_r - s / tau_d

    return np.hstack((dv, dm, dn, dh, dq, ds))
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


def tau_peak_function(tau_d, tau_r, tau_d_q):

    dt = 0.01
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
# ------------------------------------------------------------------#


def tau_d_q_function(tau_d, tau_r, tau_hat):

    # set an interval for tau_d_q
    tau_d_q_left = 1.0
    while tau_peak_function(tau_d, tau_r, tau_d_q_left) > tau_hat:
        tau_d_q_left *= 0.5

    tau_d_q_right = tau_r
    while tau_peak_function(tau_d, tau_r, tau_d_q_right) < tau_hat:
        tau_d_q_right *= 2.0

    # bisection method
    while tau_d_q_right - tau_d_q_left > 1e-12:
        tau_d_q_mid = 0.5 * (tau_d_q_left + tau_d_q_right)
        if (tau_peak_function(tau_d, tau_r, tau_d_q_mid) <= tau_hat):
            tau_d_q_left = tau_d_q_mid
        else:
            tau_d_q_right = tau_d_q_mid

    return 0.5 * (tau_d_q_left + tau_d_q_right)


tau_d_q = tau_d_q_function(tau_d, tau_r, tau_peak)
