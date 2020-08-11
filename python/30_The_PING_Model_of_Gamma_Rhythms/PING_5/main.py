from scipy.integrate import odeint
from numpy.random import randn, rand, uniform
from time import time
import numpy as np
import pylab as pl
import lib
import os

num_e = 200
num_i = 50
sigma_e = 0.05
i_ext_e = 1.4 * np.ones(num_e) * (1 + sigma_e * randn(num_e))
sigma_i = 0.00
i_ext_i = 0.0 * np.ones(num_i) * (1 + sigma_i * randn(num_i))
g_hat_ee = 0
g_hat_ei = 0.25
g_hat_ie = 0.25
g_hat_ii = 0.25
p_ee = 0.5
p_ei = 0.5
p_ie = 0.5
p_ii = 0.5

v_rev_e = 0.0
v_rev_i = -75.0
tau_r_e = 0.5
tau_peak_e = 0.5
tau_d_e = 3.0
tau_r_i = 0.5
tau_peak_i = 0.5
tau_d_i = 9.0
t_final = 200.0
dt = 0.02

spikeThreshold = -20
# Process network parameters a bit:

u_ee = rand(num_e, num_e)
u_ei = rand(num_e, num_i)
u_ie = rand(num_i, num_e)
u_ii = rand(num_i, num_i)
g_ee = (g_hat_ee * (u_ee < p_ee) / (num_e * p_ee)).T
g_ei = (g_hat_ei * (u_ei < p_ei) / (num_e * p_ei)).T
g_ie = (g_hat_ie * (u_ie < p_ie) / (num_i * p_ie)).T
g_ii = (g_hat_ii * (u_ii < p_ii) / (num_i * p_ii)).T

"""
Consider, for example, the i-th e-cell and the j-th i-cell. the
probability that there is a synaptic connection at all from the
i-th e-cell to the j-th i-cell is p_ei. if there is such a
connection, its strength is g_hat_ei/(num_e*p_ei). Note that
num_e*p_ei is the expected number of excitatory inputs into an
inhibitory cell. Therefore dividing by this quantity has the
effect that the expected value of the total excitatory
conductance affecting an inhibitory cell is g_hat_ei. 
"""

# ------------------------------------------------------------------#


if __name__ == "__main__":

    start = time()

    # os.setenv('OMP_NUM_THREADS', '1')
    # os.environ["OMP_NUM_THREADS"] = "1"

    initialVec = lib.splayState(i_ext_e,
                                rand(num_e),
                                lib.derivative)

    v_e = initialVec[:, 0]
    m_e = lib.m_e_inf(v_e)
    h_e = initialVec[:, 1]
    n_e = initialVec[:, 2]
    q_e = np.zeros(num_e)
    s_e = np.zeros(num_e)

    v_i = -75.0 * np.ones(num_i)
    m_i = lib.m_i_inf(v_i)
    h_i = lib.h_i_inf(v_i)
    n_i = lib.n_i_inf(v_i)
    q_i = np.zeros(num_i)
    s_i = np.zeros(num_i)

    initialConditions = np.hstack((v_e, h_e, n_e, q_e, s_e,
                                   v_i, h_i, n_i, q_i, s_i))
    t = np.arange(0, t_final, dt)

    sigma_e = 0.0
    sol = odeint(lib.derivativePopulation,
                 initialConditions,
                 t)

    lfp = np.mean(sol[:, :num_e], axis=1)

    t_e_spikes = []
    t_i_spikes = []
    for i in range(num_e):
        ts_e = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_e_spikes.append(ts_e)
    index = 5 * num_e
    for i in range(index, index + num_i):
        ts_i = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_i_spikes.append(ts_i)

    lib.display_time(time() - start)

    lib.spikeToFile(t_e_spikes, "t_e_spikes1.txt")
    lib.spikeToFile(t_i_spikes, "t_i_spikes1.txt")
    np.savetxt("lfp1.txt", zip(t, lfp), fmt="%18.6f")

    # --------------------------------------------------------------#
    p_ee = 1.0 / num_e
    p_ei = 1.0 / num_i
    p_ie = 1.0 / num_i
    p_ii = 1.0 / num_i

    u_ee = rand(num_e, num_e)
    u_ei = rand(num_e, num_i)
    u_ie = rand(num_i, num_e)
    u_ii = rand(num_i, num_i)
    g_ee = (g_hat_ee * (u_ee < p_ee) / (num_e * p_ee)).T
    g_ei = (g_hat_ei * (u_ei < p_ei) / (num_e * p_ei)).T
    g_ie = (g_hat_ie * (u_ie < p_ie) / (num_i * p_ie)).T
    g_ii = (g_hat_ii * (u_ii < p_ii) / (num_i * p_ii)).T

    sol = odeint(lib.derivativePopulation,
                 initialConditions,
                 t)

    lfp = np.mean(sol[:, :num_e], axis=1)

    t_e_spikes = []
    t_i_spikes = []
    for i in range(num_e):
        ts_e = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_e_spikes.append(ts_e)
    index = 5 * num_e
    for i in range(index, index + num_i):
        ts_i = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_i_spikes.append(ts_i)

    lib.display_time(time() - start)

    lib.spikeToFile(t_e_spikes, "t_e_spikes2.txt")
    lib.spikeToFile(t_i_spikes, "t_i_spikes2.txt")
    np.savetxt("lfp2.txt", zip(t, lfp), fmt="%18.6f")

    # --------------------------------------------------------------#

    # Define new network, with each cell having a fixed number of pre-synaptic
    # E-cells and I-cells:

    ni = 1  # number of excitatory and inhibitory inputs per cell

    g_ee = np.zeros((num_e, num_e))
    g_ei = np.zeros((num_e, num_i))
    g_ie = zeros((num_i, num_e))
    g_ii = np.zeros((num_i, num_i))

    i_vec = np.zeros(num_e)
    
    for j in range(num_i):
        i_vec[0] = ceil(rand() * num_e)
        for m in range(2, ni-1, -1):

            i0 = i_vec[0]
            while max(abs(i_vec - i0) == 0):
                i0 = ceil(rand()*num_e)
            i_vec[m] = i0

        g_ei[i_vec, j] = g_hat_ei / (num_e * p_ei)

    # for i=1:num_e
    for i in range(num_e):
        j_vec[0] = ceil(rand() * num_i)

        for m in range(2, ni - 1, -1):
            j0 = j_vec[0]

            while max(abs(j_vec - j0) == 0):

                j0 = ceil(rand() * num_i)

            j_vec[m] = j0

        g_ie[j_vec, i] = g_hat_ie / (num_i * p_ie)

    for j in range(num_i):  # j=1:num_i
        j_vec[0] = ceil(rand() * num_i)

        for m in range(2, ni-1, -1):
            j0 = j_vec[0]

            while max(abs(j_vec-j0) == 0):
                j0 = ceil(rand() * num_i)
            j_vec[m] = j0

        g_ii[j_vec, j] = g_hat_ii / (num_i * p_ii)



    sol = odeint(lib.derivativePopulation,
                 initialConditions,
                 t)

    lfp = np.mean(sol[:, :num_e], axis=1)

    t_e_spikes = []
    t_i_spikes = []
    for i in range(num_e):
        ts_e = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_e_spikes.append(ts_e)
    index = 5 * num_e
    for i in range(index, index + num_i):
        ts_i = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_i_spikes.append(ts_i)

    lib.display_time(time() - start)

    lib.spikeToFile(t_e_spikes, "t_e_spikes3.txt")
    lib.spikeToFile(t_i_spikes, "t_i_spikes3.txt")
    np.savetxt("lfp3.txt", zip(t, lfp), fmt="%18.6f")