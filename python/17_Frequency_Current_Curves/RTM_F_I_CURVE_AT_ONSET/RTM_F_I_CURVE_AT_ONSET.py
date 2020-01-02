from scipy.integrate import odeint
import numpy as np
from numpy import exp
import pylab as pl
from copy import copy
from time import time


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


def derivative(x0, i_ext):
    '''
    define Traub Model
    '''
    v, m, n, h, = x0
    dv = i_ext - g_na * h * m ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dm = alpha_m(v) * (1.0 - m) - beta_m(v) * m
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return np.array([dv, dm, dn, dh])


v = -70.0
m = m_inf(v)
h = 0.7  # h_inf(v)
n = 0.6  # n_inf(v)
x0 = [v, m, n, h]


def rungeKuttaIntegrator(x, h, f, i_ext):

    k1 = h * f(x, i_ext)
    k2 = h * f(x + 0.5 * k1, i_ext)
    k3 = h * f(x + 0.5 * k2, i_ext)
    k4 = h * f(x + k3, i_ext)

    x = x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    return x


c = 1
g_k = 80
g_na = 100
g_l = 0.1
v_k = -100
v_na = 50
v_l = -67

i_ext_high = 0.1194
i_ext_low = 0.1193
i_ext_vec = np.linspace(i_ext_low, i_ext_high, 11)
t_final = 20000.0
dt = 0.05


v = -70.0
m = m_inf(v)
h = 0.7  # h_inf(v)
n = 0.6  # n_inf(v)
initialConditions = [v, m, n, h]
vThreshold = -20.0


if __name__ == "__main__":

    pl.figure(figsize=(7, 3))

    start = time()

    N = int(1000 / dt)
    num_steps = int(t_final/dt)
    frequencies = np.zeros(len(i_ext_vec))

    for direction in ["forward", "backward"]:

        if direction == "backward":
            i_ext_vec = i_ext_vec[::-1]

        for ii in range(len(i_ext_vec)):

            num_spikes = 0
            t_spikes = []
            i_ext = i_ext_vec[ii]
            v = np.zeros(num_steps)
            m = np.zeros_like(v)
            n = np.zeros_like(v)
            h = np.zeros_like(v)

            for i in range(num_steps):
                v[i], m[i], n[i], h[i] = rungeKuttaIntegrator(initialConditions,
                                                              dt,
                                                              derivative,
                                                              i_ext)

                initialConditions = copy([v[i], m[i], n[i], h[i]])
                # condition to find steady state
                if ((i % N) == 0) and (i > 0):
                    maxv = max(v[i - N:i])
                    minv = min(v[i - N:i])
                    maxm = max(m[i - N:i])
                    minm = min(m[i - N:i])
                    maxn = max(n[i - N:i])
                    minn = min(n[i - N:i])
                    maxh = max(h[i - N:i])
                    minh = min(h[i - N:i])
                    if (((maxv - minv) < 0.0001 * abs(maxv + minv)) &
                        ((maxm - minm) < 0.0001 * abs(maxm + minm)) &
                        ((maxh - minh) < 0.0001 * abs(maxh + minh)) &
                            ((maxn - minn) < 0.0001 * abs(maxn + minn))):
                        frequencies[ii] = 0.0
                        print "I =%18.6f, f =%18.6f" % (i_ext, frequencies[ii])
                        break

                # spike detection
                if (v[i-1] < vThreshold) & (v[i] >= vThreshold):
                    num_spikes += 1
                    tmp = ((i - 1) * dt * (v[i - 1] - vThreshold) +
                           i * dt * (vThreshold - v[i])) / (v[i - 1] - v[i])
                    t_spikes.append(tmp)
                    # print num_spikes

                if num_spikes == 4:
                    frequencies[ii] = 1000.0 / (t_spikes[-1] - t_spikes[-2])
                    print "I =%18.6f, f =%18.6f, t =%18.6f" % (
                        i_ext, frequencies[ii], tmp)
                    break

        # save to file
        if direction == "backward":
            np.savetxt("backward.txt", zip(
                i_ext_vec, frequencies), fmt="%20.9f")
        else:
            np.savetxt("forward.txt", zip(
                i_ext_vec, frequencies), fmt="%20.9f")

    print "Done in %10.3f" % (time() - start)

    f_forward = np.loadtxt("forward.txt")
    I = f_forward[:, 0]
    f = f_forward[:, 1]

    pl.plot(I, f, "ko", label="forward")


    # plot the solid red line
    index = np.where(f == 0)[0]
    index = max(index)
    I0, f0 = I[index], f[index]
    alpha_vec = np.arange(0, 1, 101)

    I_c_low = i_ext_vec[index]
    I_c_high = i_ext_vec[index + 1]

    C_vec = np.zeros(len(alpha_vec))
    err_vec = np.zeros(len(alpha_vec))

    for ijk in range(len(alpha_vec)):

        alpha = alpha_vec[ijk]
        I_c = I_c_low * alpha + I_c_high * (1-alpha)
        y = f0 / np.sqrt(I0 - I_c)

        # measure how constant f/sqrt(I-I_c) is, with
        # the I_c defined above.
        err_vec[ijk] = (max(y) - min(y)) / np.mean(y)
        C_vec[ijk] = np.mean(y)

    # pick the I_c that makes f/sqrt(I-I_c) as constant
    ind = np.argmin(err_vec)
    alpha = alpha_vec[ind]

    I_c = I_c_low * alpha + I_c_high * (1 - alpha)
    C = C_vec[ind]
    I = I_c + np.linspace(0, 1.0, 1001) * (i_ext_high - I_c)
    pl.plot(I, C * np.sqrt(I - I_c), '-r', lw=2)

    

    f_backward = np.loadtxt("backward.txt")
    I = f_backward[:, 0][::-1]
    f = f_backward[:, 1][::-1]

    pl.plot(I, f, "ro", fillstyle="none", markersize=8, label="backward")

    pl.xlim(min(i_ext_vec), max(i_ext_vec))
    pl.xlabel("I", fontsize=16)
    pl.ylabel("frequency", fontsize=16)
    pl.legend(fontsize=14, frameon=False)
    pl.tight_layout()
    pl.tick_params(labelsize=14)
    pl.savefig("fig_17_5.png")
