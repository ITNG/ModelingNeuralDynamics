import numpy as np
import pylab as plt
from numpy import exp
from scipy.integrate import odeint


def n_slow_inf(v): return 1./(1+exp((-20-v)/5))


def n_inf(v): return 1/(1+exp((-25-v)/5))


def m_inf(v): return 1/(1+exp((-20-v)/15))


def derivative(x0, t):
    '''
    define Model
    '''
    v, n, n_slow = x0

    dn = (n_inf(v)-n)/tau_n
    dn_slow = (n_slow_inf(v) - n_slow) / tau_n_slow
    dv = (g_na * m_inf(v) * (v_na - v) + g_k * n * (v_k - v) +
          g_k_slow * n_slow * (v_k-v) + g_l*(v_l-v) + i_ext)/cm

    return [dv, dn, dn_slow]


def initial_condition(v):

    n = n_inf(v)
    n_slow = n_slow_inf(v)

    return [v, n, n_slow]


if __name__ == "__main__":

    cm = 1
    g_na = 20
    g_k = 10
    g_l = 8
    v_na = 60
    v_k = -90
    v_l = -80
    tau_n = 0.15

    tau_n_slow = 20
    g_k_slow = 5
    i_ext = 7

    t_final = 100
    dt = 0.01
    # m_steps = round(t_final/dt)
    t = np.arange(0, t_final, dt)
    x0 = initial_condition(-70.)
    sol = odeint(derivative, x0, t)
    v = sol[:,0]


    fig, ax = plt.subplots(1, figsize=(10, 3.5))
    ax.plot(t, v, label="v", lw=2, color='k')
    ax.set_xlabel("t [ms]", fontsize=14)
    ax.set_ylabel("v [mV]", fontsize=14)
    ax.set_yticks([-100,-50, 0])
    ax.tick_params(labelsize=14)
    ax.legend(frameon=False)
    ax.margins(x=0.0)
    plt.tight_layout()
    plt.savefig("fig_19.2.png")
    plt.close()
    # plt.show()