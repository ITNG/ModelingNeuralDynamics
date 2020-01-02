from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl


def beta_h(v):
    return 5.0 / (exp(-0.1 * (v + 28.0)) + 1.0)


def beta_n(v):
    return 0.625 * exp(-(v + 44.0) / 80.0)


def beta_m(v):
    return 4.0 * exp(-(v + 60.0) / 18.0)


def alpha_h(v):
    return 0.35 * exp(-(v + 58.0) / 20.0)


def alpha_m(v):
    return 0.1 * (v + 35.0) / (1.0 - exp(-(v + 35.0) / 10.0))


def alpha_n(v):
    return 0.05 * (v + 34.0) / (1.0 - exp(-0.1 * (v + 34.0)))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))

def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))

def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))

def derivative(x0, t):

    df = np.zeros(3 * N)
    
    v = x0[:N]
    h = x0[N : 2 * N]
    n = x0[2 * N :]
    m = np.zeros(N)
    I_syn = np.zeros(N)

    for i in range(N):

        m[i] = alpha_m(v[i]) / (alpha_m(v[i]) + beta_m(v[i]))
        I_Na = g_Na * m[i] ** 3 * h[i] * (v[i] - v_na)
        I_L = g_l * (v[i] - v_l)
        I_K = g_k * n[i] ** 4 * (v[i] - v_k)

        sumj = 0.0
        for j in range(N):
            if i != j:
                I_syn[i] += (v[j] - v[i])
        
        I_syn[i] *= g_gap

        df[i] = -I_Na - I_K - I_L + I_syn[i] + i_ext[i]                        # dv
        df[i + N] = alpha_h(v[i]) * (1-h[i]) - beta_h(v[i]) * h[i]          # dh
        df[i + 2 * N] = alpha_n(v[i]) * (1 - n[i]) - beta_n(v[i]) * n[i]    # dn

    return df


# parameters
N = 2
c = 1.0
g_k = 9.0
g_Na = 35.0
g_l = 0.1
v_k = -90.0
v_na = 55.0
v_l = -65.0
i_ext = [1.0, 0.0]
t_final = 200.0
dt = 0.01
g_gap = 0.01


v = np.array([-63.0, -63.0])
# m = m_inf(v)
h = h_inf(v)
n = n_inf(v)
x0 = v.tolist() + n.tolist() + h.tolist()


if __name__ == "__main__":

    t = np.arange(0, t_final, dt)
    sol = odeint(derivative, x0, t)
    v1 = sol[:, 0]
    v2 = sol[:, 1]

    fig, ax = pl.subplots(2, figsize=(7, 5), sharex=True)
    ax[0].plot(t, v1, lw=2, c="k")
    ax[1].plot(t, v2, lw=2, c="k")


    for i in range(2):
        ax[i].set_xlim(100, max(t))
        ax[i].set_xlabel("time [ms]")
    ax[0].set_yticks(range(-100, 100, 50))
    ax[0].set_ylabel("v1 [mV]")
    ax[1].set_ylabel("v2 [mV]")
    ax[0].set_ylim(-100, 50)
    ax[1].set_ylim(-64, -62)

    pl.tight_layout()
    pl.savefig("fig_21_2.png")
    pl.close()
    # pl.show()
