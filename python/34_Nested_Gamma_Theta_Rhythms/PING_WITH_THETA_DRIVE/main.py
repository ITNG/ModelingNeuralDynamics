from scipy.integrate import odeint
from numpy.random import randn, rand, uniform
from numpy import sin, pi, ones, zeros
from time import time
import numpy as np
import pylab as pl
import lib
import os


num_e = 40
num_i = 10
sigma_e = 0.05
i_ext_e = 1.4 * ones(num_e) * (1 + sigma_e * randn(num_e))
P = 125.0
alpha = 0.8


sigma_i = 0.00
i_ext_i = 0.0 * ones(num_i) * (1 + sigma_i * randn(num_i))
g_hat_ee = 0
g_hat_ei = .25
g_hat_ie = 0.25
g_hat_ii = 0.25
p_ee = 1.0
p_ei = 1.0
p_ie = 1.0
p_ii = 1.0
# See explanation below to understand what the preceding eight
# parameters mean.

v_rev_e = 0
v_rev_i = -75
tau_r_e = 0.5
tau_peak_e = 0.5
tau_d_e = 3
tau_r_i = 0.5
tau_peak_i = 0.5
tau_d_i = 9
t_final = 1000  # Time(in ms) simulated.
dt = 0.01       # Time step used in solving the differential equations.
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
    sol = odeint(lib.derivativePopulation,
                 initialConditions,
                 t)

    lfp_v = np.mean(sol[:, :num_e], axis=1)
    lfp_s = np.mean(sol[:, 4 * num_e: 5 * num_e], axis=1) 

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

    lib.spikeToFile(t_e_spikes, "t_e_spikes.txt")
    lib.spikeToFile(t_i_spikes, "t_i_spikes.txt")
    np.savetxt("lfp_v.txt", zip(t, lfp_v), fmt="%18.6f")
    np.savetxt("lfp_s.txt", zip(t, lfp_s), fmt="%18.6f")
