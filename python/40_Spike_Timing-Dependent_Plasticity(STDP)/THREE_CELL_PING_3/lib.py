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


def derivativePopulation(x0, t, g_ee):

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


def display_time(time):
    ''' print wall time '''

    hour = int(time/3600)
    minute = (int(time % 3600))//60
    second = time-(3600.*hour+60.*minute)
    print("Done in %d hours %d minutes %.6f seconds"
          % (hour, minute, second))
# -------------------------------------------------------------------#


tau_dq_e = tau_d_q_function(tau_d_e, tau_r_e, tau_peak_e)
tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)
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


def read_from_file(fileName):

    with open(fileName, "r") as text:
        data = []
        for line in text:
            line = line.split()
            line = [float(i) for i in line]
            data.append(line)

    return data
