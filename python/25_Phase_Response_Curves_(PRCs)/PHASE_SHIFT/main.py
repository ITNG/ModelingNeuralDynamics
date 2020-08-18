from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as plt
from lib import tau_peak_function, tau_d_q_function


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


def euler_integrator(x0, f, t, dt, ts=-100):

    n = len(x0)
    nstep = len(t)
    sol = np.zeros((nstep, n))
    sol[0, :] = x0

    for i in range(1, nstep):
        x0 += f(x0, t[i], ts) * dt
        sol[i, :] = x0

    return sol


def rk4_integrator(x0, f, t, dt, ts=-100):
    n = len(x0)
    nstep = len(t)
    sol = np.zeros((nstep, n))
    sol[0, :] = x0
    for i in range(1, nstep):
        k1 = dt * f(x0, t[i], ts)
        k2 = dt * f(x0 + 0.5 * k1, t[i] + 0.5 * dt, ts)
        k3 = dt * f(x0 + 0.5 * k2, t[i] + 0.5 * dt, ts)
        k4 = dt * f(x0 + k3, t[i] + dt, ts)
        x0 += (k1 + 2.0 * (k2 + k3) + k4) / 6.0

        sol[i, :] = x0

    return sol


def ode(x0, t, ts):
    '''
    define Traub Model
    '''

    v, n, h, q, s = x0

    dv = i_ext - g_na * h * m_inf(v) ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return np.array([dv, dn, dh, 0, 0])


def ode_perturbed(x0, t, ts):
    '''
    define Traub Model
    '''

    v, n, h, q, s = x0
    if np.abs(t - ts) < (0.5*dt):
        q = 1
        print(".")
    dv = i_ext - g_na * h * m_inf(v) ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l) - \
        g_syn * s * v
    # dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h
    dq = - q / tau_d_q
    ds = q * (1 - s)/tau_r - s/tau_d

    return np.array([dv, dn, dh, dq, ds])


def initial_condition(v):
    # m = m_inf(v)
    h = h_inf(v)
    n = n_inf(v)
    q = 0.0
    s = 0.0
    return [v, n, h, q, s]


# PARAMETERS--------------------------------------------------------#
c = 1
g_k = 80
g_na = 100
g_l = 0.1
v_k = -100
v_na = 50
v_l = -67
i_ext = 0.3
t_final = 300
dt = 0.02


#  AMPA-like synaptic input pulse
tau_r = 0.5
tau_peak = 0.5
tau_d = 2.0
g_syn = 12


if __name__ == "__main__":

    fig, ax = plt.subplots(2, figsize=(7, 5), sharex=True)

    tau_d_q = tau_d_q_function(tau_d, tau_r, tau_peak, dt)
    x0 = initial_condition(-77.71)
    t = np.arange(0, t_final, dt)
    sol = rk4_integrator(x0, ode, t, dt)
    v = sol[:, 0]

    ax[0].plot(t, v, lw=3, c="b", alpha=0.5)
    ax[1].plot(t, v, lw=3, c="b", alpha=0.5)
    ax[0].set_xlim(min(t), max(t))
    ax[0].set_ylim(-100, 50)
    ax[1].set_xlabel("time [ms]")
    ax[0].set_ylabel("v [mV]")
    ax[0].set_yticks(range(-100, 100, 50))
    for i in range(2):
        ax[i].tick_params(labelsize=14)
    plt.tight_layout()

    t_stimulation = 100
    sol = rk4_integrator(x0, ode_perturbed, t, dt, ts=t_stimulation)
    v1 = sol[:, 0]
    ax[0].plot(t, v1, lw=1, c="r")
    ax[0].axvline(x=t_stimulation, ls="--", lw=1, color="k", alpha=0.5)

    t_stimulation = 120
    sol = rk4_integrator(x0, ode_perturbed, t, dt, ts=t_stimulation)
    v1 = sol[:, 0]
    ax[1].plot(t, v1, lw=1, c="r")
    ax[1].axvline(x=t_stimulation, ls="--", lw=1, color="k", alpha=0.5)

    plt.savefig("fig_5_2.png")
    plt.show()
