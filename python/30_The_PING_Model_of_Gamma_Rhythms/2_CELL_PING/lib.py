from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl
from main import *


def h_e_inf(v):
    alpha_h = 0.128 * exp(-(v + 50) / 18)
    beta_h = 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))
    return alpha_h / (alpha_h + beta_h)


def h_i_inf(v):
    alpha_h = 0.07 * exp(-(v + 58.0) / 20.0)
    beta_h = 1.0 / (exp(-0.1 * (v + 28.0)) + 1.0)
    return alpha_h / (alpha_h + beta_h)


def m_e_inf(v):
    alpha_m = 0.32 * (v + 54.0) / (1.0 - exp(-(v + 54.0) / 4.0))
    beta_m = 0.28 * (v + 27.0) / (exp((v + 27.0) / 5.0) - 1.0)
    return alpha_m / (alpha_m + beta_m)


def m_i_inf(v):
    alpha_m = 0.1 * (v + 35.0) / (1.0 - exp(-(v + 35.0) / 10.0))
    beta_m = 4.0 * exp(-(v + 60.0) / 18.0)
    return alpha_m / (alpha_m + beta_m)


def n_e_inf(v):
    alpha_n = 0.032 * (v + 52.0) / (1.0 - exp(-(v + 52.0) / 5.0))
    beta_n = 0.5 * exp(-(v + 57.0) / 40.0)
    return alpha_n / (alpha_n + beta_n)


def n_i_inf(v):
    alpha_n = -0.01 * (v + 34.0) / (exp(-0.1 * (v + 34.0)) - 1.0)
    beta_n = 0.125 * exp(-(v + 44.0) / 80.0)
    return alpha_n / (alpha_n + beta_n)


def tau_h_e(v):
    alpha_h = 0.128 * exp(-(v + 50.0) / 18.0)
    beta_h = 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))
    return 1.0 / (alpha_h + beta_h)


def tau_h_i(v):
    alpha_h = 0.07 * exp(-(v + 58.0) / 20.0)
    beta_h = 1.0 / (exp(-0.1 * (v + 28.0)) + 1.0)
    tau_h = 1.0 / (alpha_h + beta_h)
    phi = 5.0
    return tau_h / phi


def tau_n_e(v):
    alpha_n = 0.032 * (v + 52.0) / (1.0 - exp(-(v + 52.0) / 5.0))
    beta_n = 0.5 * exp(-(v + 57.0) / 40.0)
    return 1.0 / (alpha_n + beta_n)


def tau_n_i(v):
    alpha_n = -0.01 * (v + 34). / (exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * exp(-(v + 44) / 80)
    tau_n = 1. / (alpha_n + beta_n)
    phi = 5
    return tau_n / phi


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


def derivative(x0, t):

    v_e, m_e, h_e, n_e, q_e, s_e, v_i, h_i, n_i, q_i, s_i = x0

    I_Na_e = 100 * h_e * m_e_inf ** 3 * (v_e - 50.0)
    I_K_e = 80 * n_e ** 4 * (v_e + 100.0)
    I_L_e = 0.1 * (v_e + 67.0)
    I_syn_e = g_ie * s_i * (v_rev_i - v_e)

    dv_e = i_ext_e - I_L_e - I_K_e - I_Na_e + I_syn_e
    dn_e = (n_e_inf(v_e) - n_e) / tau_n_e(v_e)
    dq_e = 0.5 * (1 + np.tanh(0.1 * v_e)) * (1.0 - q_e) / 0.1 - q_e / tau_dq_e;
    
    

    dh_e = alpha_h_e(v_e) * (1.0 - h_e) - beta_h_e(v_e) * h_e
    dn_e = alpha_n_e(v_e) * (1.0 - n_e) - beta_n_e(v_e) * n_e
    dq_e = 0.5 * (1.0 + np.tanh(0.1 * v_e)) * (1.0 - q_e)/tau_r_e - q_e/tau_d_e
    ds_e = q_e * (1.0 - s_e) / tau_r_e - s_e / tau_d_e

    m_i = alpha_m_i(v_i) / (alpha_m_i(v_i) + beta_m_i(v_i))
    I_Na_i = 35.0 * m_i ** 3 * h_i * (v_i - 55.0)
    I_K_i = 9.0 * n_i ** 4 * (v_i + 90.0)
    I_L_i = 0.1 * (v_i + 65.0)
    I_syn_i = g_ei * s_e * (v_rev_e - v_i)

    dv_i = -I_Na_i - I_K_i - I_L_i + i_ext_i + I_syn_i
    dh_i = alpha_h_i(v_i) * (1.0 - h_i) - beta_h_i(v_i) * h_i
    dn_i = alpha_n_i(v_i) * (1 - n_i) - beta_n_i(v_i) * n_i
    dq_i = 0.5 * (1.0 + np.tanh(0.1 * v_i)) * \
        (1.0 - q_i) / tau_r_i - q_i / tau_d_i
    ds_i = q_i * (1.0 - s_i) / tau_r_i - s_i / tau_d_i

    return np.array([dv_e, dm_e, dh_e, dn_e, dq_e, ds_e,
                     dv_i, dh_i, dn_i, dq_i, ds_i])


tau_dq_e = tau_d_q_function(tau_d_e, tau_r_e, tau_peak_e)
tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)
