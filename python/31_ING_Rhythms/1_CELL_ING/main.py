from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl


def alpha_h(v):
    return  0.35 * exp(-(v + 58.0) / 20.0)

def alpha_m(v):
    return 0.1 * (v + 35.0) / (1.0 - exp(-0.1 * (v + 35.0)))

def alpha_n(v):
    return -0.05 * (v + 34.0) / (exp(-0.1 * (v + 34.0)) - 1.0)

def beta_h(v):
    return 5.0 / (exp(-0.1 * (v + 28.0)) + 1.0)

def beta_m(v):
    return 4.0 * exp(-(v + 60.0) / 18.0)

def beta_n(v):
    return 0.625 * exp(-(v + 44.0) / 80.0)

def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))

def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))

def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))

def derivative(x0, t):

    v, h, n, q, s = x0

    I_Na = g_Na * m_inf(v) ** 3 * h * (v - v_na)
    I_L = g_l * (v - v_l)
    I_K = g_k * n ** 4 * (v - v_k)
    I_syn = g_ii * s * (v_rev_i - v)

    dv = -I_Na - I_K - I_L + i_ext + I_syn
    dh = alpha_h(v) * (1-h) - beta_h(v) * h
    dn = alpha_n(v) * (1 - n) - beta_n(v) * n
    dq = 0.5 * (1.0 + np.tanh(0.1 * v)) * 10.0 * (1 - q) - q / tau_dq_i
    ds = q * (1.0 - s) / tau_r_i - s / tau_d_i

    return [dv, dh, dn, dq, ds]


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
        s_inc_tmp = exp(-(t + dt05) / tau_d_q) * (1.0 - s_tmp) / tau_r - s_tmp / tau_d
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


c = 1.0
g_k = 9.0
g_Na = 35.0
g_l = 0.1
v_k = -90.0
v_na = 55.0
v_l = -65.0
i_ext = 1.5
g_ii = 0.5
v_rev_i = -75.0
tau_r_i = 0.5
tau_peak_i = 0.5
tau_d_i = 9.0

t_final = 200.0
dt = 0.01


v_i = -75.0
h_i = 0.1
n_i = 0.1
q_i = 0.0
s_i = 0.0
x0 = [v_i, h_i, n_i, q_i, s_i]


if __name__ == "__main__":

    tau_dq_i = tau_d_q_function(tau_d_i, tau_r_i, tau_peak_i)


    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]

    pl.figure(figsize=(7, 3))
    pl.plot(t, v, lw=2, c="k")
    pl.xlim(min(t), max(t))
    pl.ylim(-100, 50)
    pl.xlabel("time [ms]")
    pl.ylabel("v [mV]")
    pl.yticks(range(-100, 100, 50))
    pl.tight_layout()
    pl.savefig("fig_5_3.png")
    # pl.show()
