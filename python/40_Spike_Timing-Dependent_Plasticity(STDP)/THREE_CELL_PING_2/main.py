from scipy.integrate import odeint
from numpy.random import rand
from numpy import exp
from time import time
import numpy as np
import pylab as pl
from lib import *
import lib


num_e = 2
num_i = 1
i_ext_e = [0.4, 0.8]
i_ext_i = [0]

v_rev_e = 0.0
v_rev_i = -75.0
tau_r_e = 0.5
tau_peak_e = 0.5
tau_d_e = 3.0
tau_r_i = 0.5
tau_peak_i = 0.5
tau_d_i = 9.0

t_final = 500.0
dt = 0.02

spikeThreshold = -20

g_ee = np.zeros((num_e, num_e)).T
g_ei = np.zeros((num_e, num_i)).T
g_ie = np.zeros((num_i, num_e)).T
g_ii = np.zeros((num_i, num_i)).T

g_ei[0, :] = 0.125 * np.ones(num_e)
g_ie[:, 0] = 0.25 * np.ones(num_e)
g_ii[0, 0] = 0.25


if __name__ == "__main__":

    start = time()

    v_e = -70.0 * np.ones(num_e)
    m_e = m_e_inf(v_e)
    h_e = h_e_inf(v_e)
    n_e = n_e_inf(v_e)
    q_e = np.zeros(num_e)
    s_e = np.zeros(num_e)

    v_i = -75 * np.ones(num_i)
    m_i = m_i_inf(v_i)
    h_i = h_i_inf(v_i)
    n_i = n_i_inf(v_i)
    q_i = np.zeros(num_i)
    s_i = np.zeros(num_i)

    initialConditions = np.hstack((v_e, h_e, n_e, q_e, s_e,
                                   v_i, h_i, n_i, q_i, s_i))

    t = np.arange(0, t_final, dt)
    
    g_ee[0, 1] = g_ee[1, 0] = 0.05
    sol = odeint(derivativePopulation,
                 initialConditions,
                 t,
                 args=(g_ee,))

    t_e_spikes = []
    t_i_spikes = []
    for i in range(num_e):
        ts_e = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_e_spikes.append(ts_e)
    index = 5 * num_e
    for i in range(index, index + num_i):
        ts_i = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_i_spikes.append(ts_i)

    lib.spikeToFile(t_e_spikes, "t_e_spikes1.txt")
    lib.spikeToFile(t_i_spikes, "t_i_spikes1.txt")

    g_ee[0, 1] = g_ee[1, 0] = 0.4
    sol = odeint(derivativePopulation,
                 initialConditions,
                 t,
                 args=(g_ee,))

    t_e_spikes = []
    t_i_spikes = []
    for i in range(num_e):
        ts_e = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_e_spikes.append(ts_e)
    index = 5 * num_e
    for i in range(index, index + num_i):
        ts_i = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_i_spikes.append(ts_i)

    lib.spikeToFile(t_e_spikes, "t_e_spikes2.txt")
    lib.spikeToFile(t_i_spikes, "t_i_spikes2.txt")

    lib.display_time(time() - start)
