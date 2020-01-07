from scipy.integrate import odeint
from numpy import exp, ones, zeros
from numpy.random import randn, rand
from time import time
import numpy as np
import pylab as pl
import pylab as pl
import lib


# Define strengths of h - and A-currents in O-LM cells

g_h = 12.0
g_A = 22.0

# Define network parameters:

num_e = 20
num_i = 5
num_o = 5
sigma_e = 0.05
i_ext_e = 1.8 * ones(num_e) * (1 + sigma_e * randn(num_e))
sigma_i = 0.10
i_ext_i = 1.0 * ones(num_i) * (1 + sigma_i * randn(num_i))
sigma_o = 0.05
i_ext_o = -2.0 * ones(num_o) * (1 + sigma_o * randn(num_o))

g_hat_ee = 0.00
g_hat_ei = 0.25
g_hat_eo = 0.00
g_hat_ie = 0.25
g_hat_ii = 0.25
g_hat_io = 0.50
g_hat_oe = 1.00
g_hat_oi = 0.50
g_hat_oo = 0.00

p_ee = 1.0
p_ei = 1.0
p_eo = 1.0
p_ie = 1.0
p_ii = 1.0
p_io = 1.0
p_oe = 1.0
p_oi = 1.0
p_oo = 1.0

# See explanation below to understand what the preceding
# parameters mean.

v_rev_e = 0
v_rev_i = -75.0
v_rev_o = -75.0
tau_rise_e = 0.5
tau_peak_e = 0.5
tau_d_e = 3.0
tau_rise_i = 0.5
tau_peak_i = 0.5
tau_d_i = 9.0
tau_rise_o = 0.5
tau_peak_o = 0.5
tau_d_o = 20.0

t_final = 2  # Time(in ms) simulated.
dt = 0.01  # Time step used in solving the differential equations.
spikeThreshold = -20.0

# Process network parameters a bit:

u_ee = rand(num_e, num_e)
u_ei = rand(num_e, num_i)
u_eo = rand(num_e, num_o)
u_ie = rand(num_i, num_e)
u_ii = rand(num_i, num_i)
u_io = rand(num_i, num_o)
u_oe = rand(num_o, num_e)
u_oi = rand(num_o, num_i)
u_oo = rand(num_o, num_o)

g_ee = (g_hat_ee*(u_ee < p_ee)/(num_e*p_ee)).T
g_ei = (g_hat_ei*(u_ei < p_ei)/(num_e*p_ei)).T
g_eo = (g_hat_eo*(u_eo < p_eo)/(num_e*p_eo)).T

g_ie = (g_hat_ie*(u_ie < p_ie)/(num_i*p_ie)).T
g_ii = (g_hat_ii*(u_ii < p_ii)/(num_i*p_ii)).T
g_io = (g_hat_io*(u_io < p_io)/(num_i*p_io)).T

g_oe = (g_hat_oe*(u_oe < p_oe)/(num_o*p_oe)).T
g_oi = (g_hat_oi*(u_oi < p_oi)/(num_o*p_oi)).T
g_oo = (g_hat_oo*(u_oo < p_oo)/(num_o*p_oo)).T

# Consider, for example, the i-th e-cell and the j-th i-cell. the
# probability that there is a synaptic connection at all from the
# i-th e-cell to the j-th i-cell is p_ei. if there is such a
# connection, its strength is g_hat_ei/(num_e*p_ei). Note that
# num_e*p_ei is the expected number of excitatory inputs into an
# inhibitory cell. Therefore dividing by this quantity has the
# effect that the expected value of the total excitatory
# conductance affecting an inhibitory cell is g_hat_ei.

if __name__ == "__main__":

    start = time()

    iv = lib.rtmInit(i_ext_e, rand(num_e))
    v_e = iv[:, 0]
    m_e = lib.m_e_inf(v_e)
    h_e = iv[:, 1]
    n_e = iv[:, 2]
    q_e = zeros(num_e)
    s_e = zeros(num_e)

    iv = lib.wbInit(i_ext_i, rand(num_i))
    v_i = iv[:, 0]
    m_i = lib.m_i_inf(v_i)
    h_i = iv[:, 1]
    n_i = iv[:, 2]
    q_i = zeros(num_i)
    s_i = zeros(num_i)

    iv = lib.olmInit(i_ext_o, rand(num_o))
    v_o = iv[:, 0]
    m_o = lib.m_o_inf(v_o)
    h_o = iv[:, 1]
    n_o = iv[:, 2]
    r_o = iv[:, 3]
    a_o = iv[:, 4]
    b_o = iv[:, 5]
    q_o = zeros(num_o)
    s_o = zeros(num_o)

    x0 = np.hstack((
        v_e, h_e, n_e, q_e, s_e,
        v_i, h_i, n_i, q_i, s_i,
        v_o, h_o, n_o, r_o, a_o, b_o, q_o, s_o))

    t = np.arange(0, t_final, dt)
    sol = odeint(lib.derivative, x0, t)
    lfp_v_e = np.mean(sol[:, :num_e], axis=1)
    lfp_s_e = np.mean(sol[:, 4 * num_e : 5 * num_e], axis=1)

    t_e_spikes = []
    t_i_spikes = []
    t_o_spikes = []
    for i in range(num_e):
        ts_e = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_e_spikes.append(ts_e)
    index = 5 * num_e
    for i in range(index, index + num_i):
        ts_i = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_i_spikes.append(ts_i)
    
    index = index + 5 * num_i
    for i in range(index, index + num_o):
        ts_o = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_o_spikes.append(ts_o)

    lib.display_time(time() - start)

    lib.spikeToFile(t_e_spikes, "t_e_spikes.txt")
    lib.spikeToFile(t_i_spikes, "t_i_spikes.txt")
    lib.spikeToFile(t_o_spikes, "t_o_spikes.txt")
    np.savetxt("lfp_v_e.txt", zip(t, lfp_v_e), fmt="%18.6f")
    np.savetxt("lfp_s_e.txt", zip(t, lfp_s_e), fmt="%18.6f")

    