from numpy import exp, matmul
from copy import copy
import numpy as np
from main import *


def tau_n_i(v):
    alpha_n = -0.01 * (v + 34) / (exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * exp(-(v + 44) / 80)
    tau_n = 1. / (alpha_n + beta_n)
    phi = 5.0
    tau_n /= phi
    return tau_n


def tau_h_i(v):
    alpha_h = 0.07 * exp(-(v + 58) / 20)
    beta_h = 1. / (exp(-0.1 * (v + 28)) + 1)
    tau_h = 1. / (alpha_h + beta_h)
    phi = 5.0
    tau_h = tau_h / phi
    return tau_h


def n_i_inf(v):
    alpha_n = -0.01 * (v + 34) / (exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * exp(-(v + 44) / 80)
    return alpha_n / (alpha_n + beta_n)


def m_i_inf(v):
    alpha_m = 0.1 * (v + 35) / (1 - exp(-(v + 35) / 10))
    beta_m = 4 * exp(-(v + 60) / 18)
    return alpha_m / (alpha_m + beta_m)


def h_i_inf(v):
    alpha_h = 0.07 * exp(-(v + 58) / 20)
    beta_h = 1 / (exp(-0.1 * (v + 28)) + 1)
    return alpha_h / (alpha_h + beta_h)


def derivative(x0, t):

    v = x0[:num_i]
    h = x0[num_i: 2 * num_i]
    n = x0[2 * num_i: 3 * num_i]
    q = x0[3 * num_i: 4 * num_i]
    s = x0[4 * num_i: 5 * num_i]

    I_Na = g_Na * m_i_inf(v) ** 3 * h * (v - v_na)
    I_L = g_l * (v - v_l)
    I_K = g_k * n ** 4 * (v - v_k)
    I_syn = matmul(g_ii, s) * (v_rev_i - v)  # + G_gap * v - c * v

    dv = -I_Na - I_K - I_L + i_ext + I_syn
    dh = (h_i_inf(v) - h) / tau_h_i(v)
    dn = (n_i_inf(v) - n) / tau_n_i(v)
    dq = 0.5 * (1.0 + np.tanh(0.1 * v)) * 10.0 * (1 - q) - q / tau_dq_i
    ds = q * (1.0 - s) / tau_r_i - s / tau_d_i

    return np.hstack((dv, dh, dn, dq, ds))


def derivativeSingle(x0, t=0):

    n0 = len(x0) / 3
    v, h, n = x0[:n0], x0[n0: (2 * n0)], x0[(2 * n0): (3 * n0)]

    I_Na = g_Na * m_i_inf(v) ** 3 * h * (v - v_na)
    I_L = g_l * (v - v_l)
    I_K = g_k * n ** 4 * (v - v_k)

    dv = -I_Na - I_K - I_L + i_ext
    dh = (h_i_inf(v) - h) / tau_h_i(v)
    dn = (n_i_inf(v) - n) / tau_n_i(v)

    return np.hstack((dv, dh, dn))

# -------------------------------------------------------------------#


def eulerIntegrator(x, dt, f):  

    x += f(x) * dt
    return x
# -------------------------------------------------------------------#


def splayState(i_ext, phiVec, f):
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
    N = len(i_ext)
    numSpikes = np.zeros(N, dtype=int)
    done = np.zeros(N)              # done[i]=1 indicates that we are
    # done with neuron i.
    iteration = 0
    numSteps = int(t_final / dt)

    v = -70.0 * np.ones(N)
    # m = m_i_inf(v)
    n = n_i_inf(v)
    h = h_i_inf(v)
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

        x = eulerIntegrator(x0, dt, f)
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


def display_time(time):
    ''' print wall time '''

    hour = int(time/3600)
    minute = (int(time % 3600))/60
    second = time-(3600.*hour+60.*minute)
    print "Done in %d hours %d minutes %.6f seconds" \
        % (hour, minute, second)


def spikeToFile(t_spikes, fileName):

    f = open(fileName, "w")
    n = len(t_spikes)
    for i in range(n):
        for j in range(len(t_spikes[i])):
            f.write("%18.4f" % t_spikes[i][j])
        f.write("\n")

    f.close()


def read_from_file(fileName):

    with open(fileName, "r") as text:
        data = []
        for line in text:
            line = line.split()
            line = [float(i) for i in line]
            data.append(line)

    return data


tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)
