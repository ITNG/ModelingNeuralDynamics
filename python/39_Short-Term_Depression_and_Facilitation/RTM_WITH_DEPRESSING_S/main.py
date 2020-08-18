from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl



def alpha_h(v):
    return 0.128 * exp(-(v + 50.0) / 18.0)


def alpha_m(v):
    return 0.32 * (v + 54) / (1.0 - exp(-(v + 54.0) / 4.0))


def alpha_n(v):
    return 0.032 * (v + 52) / (1.0 - exp(-(v + 52.0) / 5.0))


def beta_h(v):
    return 4.0 / (1.0 + exp(-(v + 27.0) / 5.0))


def beta_m(v):
    return 0.28 * (v + 27.0) / (exp((v + 27.0) / 5.0) - 1.0)


def beta_n(v):
    return 0.5 * exp(-(v + 57.0) / 40.0)


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def derivative(x0, t):
    '''
    define Traub Model
    '''
    v, m, n, h, p, q, s = x0
    dv = i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h
    dp = -C * (1 + np.tanh(0.1 * v)) * p * np.log(1 / (1 - U)) + (1 - p - q) / tau_rec
    dq = C * (1 + np.tanh(0.1 * v)) * p * np.log(1 / (1 - U)) - q / tau_d_q
    ds = q * (1 - s)/tau_r - s/tau_d

    return [dv, dm, dn, dh, dp, dq, ds]

def initial_condition(v):
    m = m_inf(v)
    h = h_inf(v)
    n = n_inf(v)
    p = 1.0
    q = 0.0
    s = 0.0
    return [v, m, n, h, p, q, s]


c = 1
g_k = 80
g_na = 100
g_l = 0.1
v_k = -100
v_na = 50
v_l = -67
i_ext = 0.5
t_final = 200
dt = 0.01

tau_rec = 500
tau_d_q = 5
tau_r = 3
tau_d = 9

x0 = initial_condition(-70)

U = 0.5
C = 1.45

if __name__ == "__main__":

    t = np.arange(0, t_final+dt, dt)
    sol = odeint(derivative, x0, t)
    v = sol[:, 0]
    p = sol[:, 4]
    q = sol[:, 5]
    s = sol[:, 6]

    num_spikes = 0
    for i in range(len(v) - 1):
        if (v[i] < -20) and (v[i + 1] >= -20):
            num_spikes += 1
    gamma = C * (1 + np.tanh(v / 10.0))
    integral_of_delta_function_per_period = np.sum(gamma) * dt / num_spikes
    

    fig, ax = pl.subplots(nrows=2, ncols=2, figsize=(6, 5), sharex=True)
    ax[0, 0].plot(t, v, lw=2, c="k")
    ax[0, 1].plot(t, p, lw=2, c='k')
    ax[1, 0].plot(t, q, lw=2, c='k')
    ax[1, 1].plot(t, s, lw=2, c='k')
    
    ax[0, 0].set_ylabel("v [mV]")
    ax[0, 1].set_ylabel('p')
    ax[1, 0].set_ylabel('q')
    ax[1, 1].set_ylabel('s')
    ax[1, 0].set_xlabel("time [ms]")
    ax[1, 1].set_xlabel("time [ms]")
    ax[0, 0].set_yticks(range(-100, 100, 50))
    ax[0, 0].set_ylim(-100, 50)
    ax[0, 1].set_xlim(min(t), max(t))
    ax[0, 1].set_ylim(0, 1)

    pl.tight_layout()
    
    pl.savefig("fig_39_2.png")
    # pl.show()
