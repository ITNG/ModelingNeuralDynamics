from numpy import exp, matmul
from copy import copy
import numpy as np
import pylab as pl
from main import *


def alpha_h_o(v):
    return 0.07 * exp(-(v + 63.0) / 20.0)


def alpha_m_o(v):
    q = (v+38.0)/10.0
    return q / (1.0 - exp(-q))


def alpha_n_o(v):
    return 0.018 * (v - 25.0) / (1.0 - exp(-(v - 25.0) / 25.0))


def beta_h_o(v):
    return 1.0 / (exp(-(v + 33.0) / 10.0) + 1.0)


def beta_m_o(v):
    return 4.0 * exp(-(v + 65.0) / 18.0)


def beta_n_o(v):
    return 0.0036 * (35.0 - v) / (1.0 - exp(-(35.0 - v) / 12.0))


def b_o_inf(v):
    return 1.0/(1.0+exp((v+71.0)/7.3))


def a_o_inf(v):
    return 1.0 / (1.0 + exp(-(v + 14.0) / 16.6))


def h_o_inf(v):
    return alpha_h_o(v) / (alpha_h_o(v) + beta_h_o(v))


def m_o_inf(v):
    return alpha_m_o(v) / (alpha_m_o(v) + beta_m_o(v))


def n_o_inf(v):
    return alpha_n_o(v) / (alpha_n_o(v) + beta_n_o(v))


def r_o_inf(v):
    return 1.0 / (1.0 + exp((v + 84.0) / 10.2))


def tau_a_o(v):
    return 5.0


def tau_b_o(v):
    return 1. / (0.000009 / exp((v - 26) / 28.5) +
                 0.014 / (0.2 + exp(-(v + 70.0) / 11.0)))


def tau_r_o(v):
    return 1.0 / (exp(-14.59 - 0.086 * v) + exp(-1.87 + 0.0701 * v))


def tau_h_e(v):
    alpha_h = 0.128 * exp(-(v + 50.0) / 18.0)
    beta_h = 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))
    return 1.0 / (alpha_h + beta_h)


def tau_h_i(v):
    phi = 5
    alpha_h = 0.07 * exp(-(v + 58.0) / 20.0)
    beta_h = 1.0 / (exp(-0.1 * (v + 28.0)) + 1.0)
    return 1.0 / (alpha_h + beta_h) / phi


def tau_h_o(v):
    return 1.0 / (alpha_h_o(v) + beta_h_o(v))


def tau_n_e(v):
    alpha_n = 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))
    beta_n = 0.5 * exp(-(v + 57) / 40)
    return 1.0 / (alpha_n + beta_n)


def tau_n_i(v):
    phi = 5.0
    alpha_n = -0.01 * (v + 34) / (exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * exp(-(v + 44) / 80)
    return 1.0 / (alpha_n + beta_n) / phi


def tau_n_o(v):
    return 1.0 / (alpha_n_o(v) + beta_n_o(v))


def h_e_inf(v):
    alpha_h = 0.128 * exp(-(v + 50.0) / 18.0)
    beta_h = 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))
    return alpha_h / (alpha_h + beta_h)


def h_i_inf(v):
    alpha_h = 0.07 * exp(-(v + 58.0) / 20.0)
    beta_h = 1.0 / (exp(-0.1 * (v + 28.0)) + 1.0)
    return alpha_h / (alpha_h + beta_h)


def m_e_inf(v):
    alpha_m = 0.32*(v+54.0)/(1.0-exp(-(v+54.0)/4.0))
    beta_m = 0.28*(v+27.0)/(exp((v+27.0)/5.0)-1.0)
    return alpha_m / (alpha_m + beta_m)


def n_e_inf(v):
    alpha_n = 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))
    beta_n = 0.5 * exp(-(v + 57) / 40)
    return alpha_n / (alpha_n + beta_n)


def m_i_inf(v):
    alpha_m = 0.1 * (v + 35.0) / (1 - exp(-(v + 35.0) / 10.0))
    beta_m = 4.0 * exp(-(v + 60.0) / 18.0)
    return alpha_m / (alpha_m + beta_m)


def n_i_inf(v):
    alpha_n = -0.01 * (v + 34) / (exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * exp(-(v + 44) / 80)
    return alpha_n / (alpha_n + beta_n)
# -------------------------------------------------------------------#


def derivative(x0, t):

    v_e, h_e = x0[:num_e], x0[num_e:2*num_e]
    n_e = x0[2 * num_e : 3 * num_e]
    q_e = x0[3 * num_e : 4 * num_e]
    s_e = x0[4 * num_e : 5 * num_e]
    
    n0 = 5 * num_e 
    v_i = x0[n0: n0 + num_i]
    h_i = x0[n0 + num_i : n0 + 2 * num_i]
    n_i = x0[n0 + 2 * num_i : n0 + 3 * num_i]
    q_i = x0[n0 + 3 * num_i : n0 + 4 * num_i]
    s_i = x0[n0 + 4 * num_i : n0 + 5 * num_i]
    
    n1 = n0 + 5 * num_i
    v_o = x0[n1: n1 + num_o]
    h_o = x0[n1 + num_o : n1 + 2 * num_o]
    n_o = x0[n1 + 2 * num_o : n1 + 3 * num_o]
    r_o = x0[n1 + 3 * num_o : n1 + 4 * num_o]
    a_o = x0[n1 + 4 * num_o : n1 + 5 * num_o]
    b_o = x0[n1 + 5 * num_o : n1 + 6 * num_o]
    q_o = x0[n1 + 6 * num_o : n1 + 7 * num_o]
    s_o = x0[n1 + 7 * num_o : n1 + 8 * num_o]

    I_K_e = 80.0 * n_e ** 4 * (v_e + 100.0)
    I_Na_e = 100.0 * h_e * m_i_inf(v_e) ** 3 * (v_e - 50.0)
    I_L_e = 0.1 * (v_e + 67.0)
    
    dv_e = i_ext_e - I_L_e - I_K_e - I_Na_e + \
        matmul(g_ee, s_e) * (v_e - v_rev_e) + \
        matmul(g_ie, s_i) * (v_e - v_rev_i) + \
        matmul(g_oe, s_o) * (v_e -v_rev_o)    
        
    dh_e = (h_e_inf(v_e) - h_e) / tau_h_e(v_e)
    dn_e = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)
    dq_e = 0.5 * (1.0 + np.tanh(0.1 * v_e)) * \
        10.0 * (1.0 - q_e) - q_e / tau_dq_e
    ds_e = q_e * (1.0 - s_e) / tau_rise_e - s_e / tau_d_e

    
    I_K_i = 9.0 * n_i ** 4 * (v_i + 90.0)
    I_Na_i = 35.0 * m_i_inf(v_i) ** 3 * h_i * (v_i - 55.0)
    I_L_i = 0.1 * (v_i + 65.0)

    dv_i = i_ext_i - I_L_i - I_K_i - I_Na_i + \
        matmul(g_ei, s_e) * (v_rev_e - v_i) + \
        matmul(g_ii, s_i) * (v_rev_i - v_i) + \
        matmul(g_oi, s_o) * (v_rev_o - v_i)
    dh_i = (h_i_inf(v_i) - h_i) / tau_h_i(v_i)
    dn_i = (n_i_inf(v_i) - n_i) / tau_n_i(v_i)
    dq_i = 0.5 * (1.0 + np.tanh(0.1 * v_i)) * \
        10.0 * (1.0 - q_i) - q_i / tau_dq_i
    ds_i = q_i * (1.0 - s_i) / tau_rise_i - s_i / tau_d_i

    I_K_o = 23.0 * n_o ** 4 * (v_o + 100.0)
    I_Na_o = 30.0 * h_o * m_o_inf(v_o) ** 3 * (v_o - 90.0)
    I_L_o = 0.05 * (v_o +70.0)
    I_H_o = g_h * r_o * (v_o + 32.9)
    I_A_o = g_A * a_o * b_o * (v_o + 90.0)

    
    dv_o = (i_ext_o - I_L_o - I_K_o - I_Na_o - I_H_o - I_A_o + \
        matmul(g_eo, s_e) * (v_rev_e - v_o) + \
        matmul(g_io, s_i) * (v_rev_i - v_o) + \
        matmul(g_oo, s_o) * (v_rev_o - v_o)) / 1.3
    
    dh_o = (h_o_inf(v_o) - h_o) / tau_h_o(v_o)
    dn_o = (n_o_inf(v_o) - n_o) / tau_n_o(v_o)
    dr_o = (r_o_inf(v_o) - r_o) / tau_r_o(v_o)
    da_o = (a_o_inf(v_o) - a_o) / tau_a_o(v_o)
    db_o = (b_o_inf(v_o) - b_o) / tau_b_o(v_o)
    dq_o = 0.5 * (1 + np.tanh(0.1 * v_o)) * (1 - q_o) / 0.1 - \
         q_o / tau_dq_o
    ds_o = q_o * (1.0 - s_o) / tau_rise_o - s_o / tau_d_o

    return np.hstack((dv_e, dh_e, dn_e, dq_e, ds_e,
            dv_i, dh_i, dn_i, dq_i, ds_i,
            dv_o, dh_o, dn_o, dr_o, da_o, db_o, dq_o, ds_o))

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


def rungeKuttaIntegrator(x, dt, f):

    k1 = dt * f(x)
    k2 = dt * f(x + 0.5 * k1)
    k3 = dt * f(x + 0.5 * k2)
    k4 = dt * f(x + k3)

    x = x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    return x
# -------------------------------------------------------------------#


def rtmInit(i_ext, phiVec):
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

    def derivativeRTM(x, t=0):

        n0 = len(x0) / 3
        v_e, h_e, n_e = x0[:n0], x0[n0: (2 * n0)], x0[(2 * n0): (3 * n0)]

        I_L_e = g_l * (v_e - v_l)
        I_K_e = g_k * n_e ** 4 * (v_e - v_k)
        I_Na_e = g_na * h_e * m_e_inf(v_e) ** 3 * (v_e - v_na)

        dv_e = i_ext - I_L_e - I_K_e - I_Na_e
        dh_e = (h_e_inf(v_e) - h_e) / tau_h_e(v_e)
        dn_e = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)

        return np.hstack((dv_e, dh_e, dn_e))

    c = 1.0
    g_k = 80.0
    g_na = 100.0
    g_l = 0.1
    v_k = -100.0
    v_na = 50.0
    v_l = -67.0

    t = 0.0
    dt = 0.01
    t_final = 2000.0
    # if fewer than max_spikes spikes occur by this time, the program gives
    # up and sets (v,h,n) equal to the values at time t_final.
    maxNumSpikes = 3
    N = len(i_ext)
    numSpikes = np.zeros(N, dtype=int)
    done = np.zeros(N)              # done[i]=1 indicates that we are
                                    # done with neuron i.
    numSteps = int(t_final / dt)

    v = -70.0 * np.ones(N)
    m = m_e_inf(v)
    n = n_e_inf(v)
    h = h_e_inf(v)
    x0 = np.hstack((v, h, n))

    tSpikes = np.zeros((N, maxNumSpikes))
    initialConition = np.zeros((N, 3))

    i = 1
    while (np.sum(done) < N) and (t < t_final):

        v_old, h_old, n_old, t_old = v, h, n, t
        x = rungeKuttaIntegrator(x0, dt, derivativeRTM)
        i += 1

        v, h, n = x[:N], x[N: (2 * N)], x[(2 * N):]
        t = i * dt
        x0 = copy(x)

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

    return initialConition
# -------------------------------------------------------------------#


def wbInit(i_ext, phiVec):
    """
    input: i_ext=column vector of external drives 
           phi_vec = column vector of phases at which 
                     neurons are to be initialized
           The length, num, of i_ext is the total number 
           of neurons.
    output: a num-by-3 array called wb_init. The columns
            contain values of v, h, and n. If i_ext(i) is below
            the firing threshold, then the i-th row of the
            matrix wb_init contains the stable equilibrium point.
            If i_ext(i) is above the firing threshold, then the
            i-th row of wb_init is a point (v,h,n) on the limit
            cycle, at phase phi_vec(i). 
    """

    def derivativeWB(x, t=0):

        n0 = len(x0) / 3
        v_i, h_i, n_i = x0[:n0], x0[n0: (2 * n0)], x0[(2 * n0): (3 * n0)]

        I_L_i = g_l * (v_i - v_l)
        I_K_i = g_k * n_i ** 4 * (v_i - v_k)
        I_Na_i = g_na * h_i * m_i_inf(v_i) ** 3 * (v_i - v_na)

        dv_i = i_ext - I_L_i - I_K_i - I_Na_i
        dh_i = (h_i_inf(v_i) - h_i) / tau_h_i(v_i)
        dn_i = (n_i_inf(v_i) - n_i) / tau_n_i(v_i)

        return np.hstack((dv_i, dh_i, dn_i))


    c = 1
    g_k = 9
    g_na = 35
    g_l = 0.1
    v_k = -90
    v_na = 55
    v_l = -65
    t = 0.0
    dt = 0.01
    t_final = 2000.0
    # if fewer than max_spikes spikes occur by this time, the program gives
    # up and sets (v,h,n) equal to the values at time t_final.

    maxNumSpikes = 3
    N = len(i_ext)
    numSpikes = np.zeros(N, dtype=int)
    done = np.zeros(N)              # done[i]=1 indicates that we are
                                    # done with neuron i.
    numSteps = int(t_final / dt)

    v = -70.0 * np.ones(N)
    m = m_i_inf(v)
    n = n_i_inf(v)
    h = h_i_inf(v)
    x0 = np.hstack((v, h, n))

    tSpikes = np.zeros((N, maxNumSpikes))
    initialConition = np.zeros((N, 3))

    i = 1
    while (np.sum(done) < N) and (t < t_final):

        v_old, h_old, n_old, t_old = v, h, n, t
        x = rungeKuttaIntegrator(x0, dt, derivativeWB)
        i += 1

        v, h, n = x[:N], x[N: (2 * N)], x[(2 * N):]
        t = i * dt
        x0 = copy(x)

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

    return initialConition

# -------------------------------------------------------------------#
def olmInit(i_ext_o, phiVec):
    """
    input: i_ext=column vector of external drives 
           phi_vec = column vector of phases at which 
                     neurons are to be initialized
           The length, num, of i_ext is the total number 
           of neurons.
    output: a num-by-6 array called olm_init. The columns
            contain values of v, h, n, r, a, b. If i_ext(i) is below
            the firing threshold, then the i-th row of the
            matrix olm_init contains the stable equilibrium point.
            If i_ext(i) is above the firing threshold, then the
            i-th row of olm_init is a point (v,h,n,r,a,b) on the limit
            cycle, at phase phi_vec(i). 
    """
    def derivativeOLM(x0, t=0):

        v_o = x0[:num_o]
        h_o = x0[num_o: 2 * num_o]
        n_o = x0[2 * num_o : 3 * num_o]
        r_o = x0[3 * num_o : 4 * num_o]
        a_o = x0[4 * num_o : 5 * num_o]
        b_o = x0[5 * num_o : 6 * num_o]

        I_K_o = g_k * n_o ** 4 * (v_o - v_k)
        I_Na_o = g_na * h_o * m_o_inf(v_o) ** 3 * (v_o - v_na)
        I_L_o = g_l * (v_o - v_l)
        I_H_o = g_h * r_o * (v_o - v_h)
        I_A_o = g_A * a_o * b_o * (v_o - v_A)
        
        dv_o = (i_ext_o - I_L_o - I_K_o - I_Na_o - I_H_o - I_A_o) / 1.3
        
        dh_o = (h_o_inf(v_o) - h_o) / tau_h_o(v_o)
        dn_o = (n_o_inf(v_o) - n_o) / tau_n_o(v_o)
        dr_o = (r_o_inf(v_o) - r_o) / tau_r_o(v_o)
        da_o = (a_o_inf(v_o) - a_o) / tau_a_o(v_o)
        db_o = (b_o_inf(v_o) - b_o) / tau_b_o(v_o)


        return np.hstack((dv_o, dh_o, dn_o, dr_o, da_o, db_o))


    c = 1.3
    g_k = 23
    g_na = 30
    g_l = 0.05
    v_k = -100
    v_na = 90
    v_l = -70
    g_h = 12
    g_A = 22
    v_h = -32.9
    v_A = -90

    max_spikes = 3
    t_final = 2000.0;
    dt = 0.01
    t = 0.0
    # if fewer than max_spikes spikes occur
    # by this time, the program gives
    # up and sets (v,h,n,r,a,b) equal 
    # to the values at time t_final.

    maxNumSpikes = 3
    N = len(i_ext_o)
    numSpikes = np.zeros(N, dtype=int)
    done = np.zeros(N)              # done[i]=1 indicates that we are
                                    # done with neuron i.
    numSteps = int(t_final / dt)

    v = -70.0 * np.ones(N)
    h = h_o_inf(v)
    n = n_o_inf(v)
    r = r_o_inf(v)
    a = a_o_inf(v)
    b = b_o_inf(v)
    x0 = np.hstack((v, h, n, r, a, b))

    tSpikes = np.zeros((N, maxNumSpikes))
    initialConition = np.zeros((N, 6))

    xx = []

    i = 1
    while (np.sum(done) < N) and (t < t_final):

        v_old = v
        h_old = h
        n_old = n
        r_old = r
        a_old = a
        b_old = b
        t_old = t
        x = rungeKuttaIntegrator(x0, dt, derivativeOLM)
        i += 1

        v = x[:N]
        h = x[N: 2 * N]
        n = x[2 * N: 3 * N]
        r = x[3 * N: 4 * N]
        a = x[4 * N: 5 * N]
        b = x[5 * N: 6 * N]
        t = i * dt
        x0 = copy(x)

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
            initialConition[k, 3] = (
                r_old[k] * (t - thr[k]) + r[k] * (thr[k] - t_old)) / dt
            initialConition[k, 4] = (
                a_old[k] * (t - thr[k]) + a[k] * (thr[k] - t_old)) / dt
            initialConition[k, 5] = (
                b_old[k] * (t - thr[k]) + b[k] * (thr[k] - t_old)) / dt 
        done[indices] = 1

    indices = np.where(done == 0)[0]
    initialConition[indices, 0] = v[indices]
    initialConition[indices, 1] = h[indices]
    initialConition[indices, 2] = n[indices]
    initialConition[indices, 3] = r[indices]
    initialConition[indices, 4] = a[indices]
    initialConition[indices, 5] = b[indices]

    return initialConition

# -------------------------------------------------------------------#

def read_from_file(fileName):

    with open(fileName, "r") as text:
        data = []
        for line in text:
            line = line.split()
            line = [float(i) for i in line]
            data.append(line)

    return data
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


def display_time(time):
    ''' print wall time '''

    hour = int(time/3600)
    minute = (int(time % 3600))/60
    second = time-(3600.*hour+60.*minute)
    print "Done in %d hours %d minutes %.6f seconds" \
        % (hour, minute, second)
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


tau_dq_e = tau_d_q_function(tau_d_e, tau_rise_e, tau_peak_e)
tau_dq_i = tau_d_q_function(tau_d_i, tau_rise_i, tau_peak_i)
tau_dq_o = tau_d_q_function(tau_d_o, tau_rise_o, tau_peak_o)
