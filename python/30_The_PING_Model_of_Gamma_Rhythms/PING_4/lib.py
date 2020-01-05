from scipy.integrate import odeint
from numpy import exp, matmul
from copy import copy
import numpy as np
import pylab as pl
from main import *


def h_e_inf(v):
    alpha_h = 0.128 * exp(-(v + 50) / 18)
    beta_h = 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))
    return (alpha_h / (alpha_h + beta_h))


def h_i_inf(v):
    alpha_h = 0.07 * exp(-(v + 58.0) / 20.0)
    beta_h = 1.0 / (exp(-0.1 * (v + 28.0)) + 1.0)
    return (alpha_h / (alpha_h + beta_h))


def m_e_inf(v):
    alpha_m = 0.32 * (v + 54.0) / (1.0 - exp(-(v + 54.0) / 4.0))
    beta_m = 0.28 * (v + 27.0) / (exp((v + 27.0) / 5.0) - 1.0)
    return (alpha_m / (alpha_m + beta_m))


def m_i_inf(v):
    alpha_m = 0.1 * (v + 35.0) / (1.0 - exp(-(v + 35.0) / 10.0))
    beta_m = 4.0 * exp(-(v + 60.0) / 18.0)
    return (alpha_m / (alpha_m + beta_m))


def n_e_inf(v):
    alpha_n = 0.032 * (v + 52.0) / (1.0 - exp(-(v + 52.0) / 5.0))
    beta_n = 0.5 * exp(-(v + 57.0) / 40.0)
    return (alpha_n / (alpha_n + beta_n))


def n_i_inf(v):
    alpha_n = -0.01 * (v + 34.0) / (exp(-0.1 * (v + 34.0)) - 1.0)
    beta_n = 0.125 * exp(-(v + 44.0) / 80.0)
    return (alpha_n / (alpha_n + beta_n))


def tau_h_e(v):
    alpha_h = 0.128 * exp(-(v + 50.0) / 18.0)
    beta_h = 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))
    return (1.0 / (alpha_h + beta_h))


def tau_h_i(v):
    alpha_h = 0.07 * exp(-(v + 58.0) / 20.0)
    beta_h = 1.0 / (exp(-0.1 * (v + 28.0)) + 1.0)
    tau_h = 1.0 / (alpha_h + beta_h)
    phi = 5.0
    return (tau_h / phi)


def tau_n_e(v):
    alpha_n = 0.032 * (v + 52.0) / (1.0 - exp(-(v + 52.0) / 5.0))
    beta_n = 0.5 * exp(-(v + 57.0) / 40.0)
    return (1.0 / (alpha_n + beta_n))


def tau_n_i(v):
    alpha_n = -0.01 * (v + 34.0) / (exp(-0.1 * (v + 34.0)) - 1.0)
    beta_n = 0.125 * exp(-(v + 44.0) / 80.0)
    tau_n = 1.0 / (alpha_n + beta_n)
    phi = 5.0
    return (tau_n / phi)
# -------------------------------------------------------------------#


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
# -------------------------------------------------------------------#


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
# -------------------------------------------------------------------#


def derivativePopulation(x0, t=0):

    v_e = x0[:num_e]
    h_e = x0[num_e: 2 * num_e]
    n_e = x0[2 * num_e: 3 * num_e]
    q_e = x0[3 * num_e: 4 * num_e]
    s_e = x0[4 * num_e: 5 * num_e]
    n = 5 * num_e
    v_i = x0[n: n + num_i]
    h_i = x0[n + num_i: n + 2 * num_i]
    n_i = x0[n + 2 * num_i: n + 3 * num_i]
    q_i = x0[n + 3 * num_i: n + 4 * num_i]
    s_i = x0[n + 4 * num_i:]

    I_L_e = 0.1 * (v_e + 67.0)
    I_K_e = 80 * n_e ** 4 * (v_e + 100.0)
    I_Na_e = 100 * h_e * m_e_inf(v_e) ** 3 * (v_e - 50.0)
    I_syn_e = matmul(g_ee, s_e) * (v_rev_e - v_e) + \
        matmul(g_ie, s_i) * (v_rev_i - v_e)

    dv_e = i_ext_e - I_L_e - I_K_e - I_Na_e + I_syn_e
    dh_e = (h_e_inf(v_e) - h_e) / tau_h_e(v_e)
    dn_e = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)
    dq_e = 0.5 * (1 + np.tanh(0.1 * v_e)) * (1.0 - q_e) * 10.0 - q_e / tau_dq_e
    ds_e = q_e * (1.0 - s_e) / tau_r_e - s_e / tau_d_e

    I_L_i = 0.1 * (v_i + 65.0)
    I_K_i = 9.0 * n_i ** 4 * (v_i + 90.0)
    I_Na_i = 35.0 * m_i_inf(v_i) ** 3 * h_i * (v_i - 55.0)
    I_syn_i = matmul(g_ei, s_e) * (v_rev_e - v_i) + \
        matmul(g_ii, s_i) * (v_rev_i - v_i)

    dv_i = i_ext_i - I_Na_i - I_K_i - I_L_i + I_syn_i
    dh_i = (h_i_inf(v_i) - h_i) / tau_h_i(v_i)
    dn_i = (n_i_inf(v_i) - n_i) / tau_n_i(v_i)
    dq_i = 0.5 * (1.0 + np.tanh(0.1 * v_i)) * (1.0 - q_i) * 10 - q_i / tau_dq_i
    ds_i = q_i * (1.0 - s_i) / tau_r_i - s_i / tau_d_i

    return np.hstack((dv_e, dh_e, dn_e, dq_e, ds_e,
                      dv_i, dh_i, dn_i, dq_i, ds_i))
# -------------------------------------------------------------------#


def derivative(x0, t=0):

    n0 = len(x0) / 3
    v_e, h_e, n_e = x0[:n0], x0[n0: (2 * n0)], x0[(2 * n0): (3 * n0)]

    I_L_e = 0.1 * (v_e + 67.0)
    I_K_e = 80 * n_e ** 4 * (v_e + 100.0)
    I_Na_e = 100 * h_e * m_e_inf(v_e) ** 3 * (v_e - 50.0)

    dv_e = i_ext_e - I_L_e - I_K_e - I_Na_e
    dh_e = (h_e_inf(v_e) - h_e) / tau_h_e(v_e)
    dn_e = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)

    return np.hstack((dv_e, dh_e, dn_e))
# -------------------------------------------------------------------#


def splayState(i_ext_e, phiVec, f):
    """
    input: i_ext=column vector of external drives 
           phi_vec = column vector of phases at which 
                     neurons are to be initialized
           The length, num, of i_ext is the total number 
           of neurons.
    output: a num-by-3 array called rtm_init. The columns
            contain values of v, h, and n. If i_ext(i) is below
            the firing threshold, then the i-th row of the
            matrix rtm_init contains the stable equilibrium point.
            If i_ext(i) is above the firing threshold, then the
            i-th row of rtm_init is a point (v,h,n) on the limit
            cycle, at phase phi_vec(i).
    """

    t = 0.0
    dt = 0.01
    t_final = 2000.0                # if fewer than max_spikes spikes occur
    # by this time, the program gives
    # up and sets (v,h,n) equal
    # to the values at time t_final.
    maxNumSpikes = 3
    N = len(i_ext_e)
    numSpikes = np.zeros(N, dtype=int)
    done = np.zeros(N)              # done[i]=1 indicates that we are
    # done with neuron i.
    iteration = 0
    numSteps = int(t_final / dt)

    v = -70.0 * np.ones(N)
    m = m_e_inf(v)
    n = n_e_inf(v)
    h = h_e_inf(v)
    x0 = np.hstack((v, h, n))

    tSpikes = np.zeros((N, maxNumSpikes))
    initialConition = np.zeros((N, 3))

    # ofile = open("v.txt", "w")
    i = 1
    while (np.sum(done) < N) and (t < t_final):

        v_old = v
        h_old = h
        n_old = n
        t_old = t

        x = rungeKuttaIntegrator(x0, dt, f)
        i += 1

        v = x[:N]
        h = x[N: (2 * N)]
        n = x[(2 * N):]
        t = i * dt

        x0 = copy(x)

        # ofile.write("%18.4f %18.4f \n" % (i*dt, v))
        indices = np.where((v_old >= spikeThreshold) &
                           (v < spikeThreshold))[0]

        numInstantSpikes = len(indices)
        for k in indices:
            numSpikes[k] += 1
            ts = (t_old * (v_old[k] - spikeThreshold) +
                  t * (spikeThreshold - v[k])) / (v_old[k] - v[k])

            if numSpikes[k] < 4:
                tSpikes[k, numSpikes[k] - 1] = ts

        thr = tSpikes[:, -1] + phiVec * (tSpikes[:, -1] - tSpikes[:, -2])
        indices = np.where((numSpikes == maxNumSpikes) &
                           (t > thr) &
                           (t_old <= thr))

        for i0 in range(len(indices)):
            k = indices[i0]
            initialConition[k, 0] = (
                v_old[k] * (t - thr[k]) + v[k] * (thr[k] - t_old)) / dt
            initialConition[k, 1] = (
                h_old[k] * (t - thr[k]) + h[k] * (thr[k] - t_old)) / dt
            initialConition[k, 2] = (
                n_old[k] * (t - thr[k]) + n[k] * (thr[k] - t_old)) / dt
        done[indices] = 1

    indices = np.where(done == 0)[0]
    initialConition[indices, 0] = v[indices]
    initialConition[indices, 1] = h[indices]
    initialConition[indices, 2] = n[indices]

    # print ("Period is %10.3f [ms]" % T)

    return initialConition
# -------------------------------------------------------------------#


def spikeDetection(t, V, spikeThreshold):
    tSpikes = []
    v = np.asarray(V)
    nSteps = len(V)

    for i in range(1, nSteps):
        if (V[i - 1] <= spikeThreshold) & (V[i] > spikeThreshold):

            ts = ((i - 1) * dt * (V[i - 1] - spikeThreshold) +
                  i * dt * (spikeThreshold - V[i])) / (V[i - 1] - V[i])
            tSpikes.append(ts)
    return tSpikes
# -------------------------------------------------------------------#


tau_dq_e = tau_d_q_function(tau_d_e, tau_r_e, tau_peak_e)
tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)
# -------------------------------------------------------------------#


def display_time(time):
    ''' print wall time '''

    hour = int(time/3600)
    minute = (int(time % 3600))/60
    second = time-(3600.*hour+60.*minute)
    print "Done in %d hours %d minutes %.6f seconds" \
        % (hour, minute, second)
# -------------------------------------------------------------------#


def spikeToFile(t_spikes, fileName):

    f = open(fileName, "w")
    n = len(t_spikes)
    for i in range(n):
        for j in range(len(t_spikes[i])):
            f.write("%18.4f" % t_spikes[i][j])
        f.write("\n")

    f.close()
# -------------------------------------------------------------------#


def rungeKuttaIntegrator(x, dt, f):

    k1 = dt * f(x)
    k2 = dt * f(x + 0.5 * k1)
    k3 = dt * f(x + 0.5 * k2)
    k4 = dt * f(x + k3)

    x = x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    return x
# -------------------------------------------------------------------#


def read_from_file(fileName):

    with open(fileName, "r") as text:
        data = []
        for line in text:
            line = line.split()
            line = [float(i) for i in line]
            data.append(line)

    return data
