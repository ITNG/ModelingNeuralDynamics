from numpy.random import rand, randn
from scipy.integrate import odeint
import networkx as nx
from time import time
import numpy as np
import lib

seed = 124875

c = 1.0
g_k = 9.0
g_Na = 35.0
g_l = 0.1
v_k = -90.0
v_na = 55.0
v_l = -65.0
i_ext = 1.5
v_rev_i = -75.0
tau_r_i = 0.5
tau_peak_i = 0.5
tau_d_i = 9.0

# Define network parameters:

num_i = 100
sigma_i = 0.0
i_ext_i = 1.5 * (1 + randn(num_i) * sigma_i)
g_hat_ii = 0.5
p_ii = 1
g_hat_gap = 0.0
p_gap = 1
# Individual gap junctions are of strength
# g_hat_gap/(p_gap*(num_i-1)).
# Any two different I-cells are gap-junctionally
# connected with probability p_gap.

t_final = 500.0   # Time (in ms) simulated.
dt = 0.01       # Time step used in solving the differential equations.
spikeThreshold = -20.0
# Process network parameters:

u_ii = rand(num_i, num_i)
g_ii = (g_hat_ii*(u_ii < p_ii)/(num_i*p_ii)).T

if __name__ == "__main__":

    start = time()

    adjacency_matrix = nx.to_numpy_array(
        nx.gnp_random_graph(num_i, p_gap, seed=seed))

    # G_gap = adjacency_matrix * g_hat_gap / (p_gap * (num_i - 1.0))
    # c = np.sum(G_gap, axis=1)

    initialConditions = lib.splayState(i_ext_i, num_i, lib.derivativeSingle)
    # v = initialConditions[:, 0]
    v = np.random.uniform(-100, 50, num_i)
    n = np.random.rand(num_i)
    h = np.random.rand(num_i)
    
    m = lib.m_i_inf(v)
    # h = initialConditions[:, 1]
    # n = initialConditions[:, 2]
    q = np.zeros(num_i)
    s = np.zeros(num_i)

    x0 = np.hstack((v, h, n, q, s))
    t = np.arange(0, t_final, dt)
    sol = odeint(lib.derivative, x0, t)

    lfp = np.mean(sol[:, :num_i], axis=1)

    t_i_spikes = []
    for i in range(num_i):
        ts_i = lib.spikeDetection(t, sol[:, i], spikeThreshold)
        t_i_spikes.append(ts_i)


    lib.display_time(time() - start)
    lib.spikeToFile(t_i_spikes, "t_i_spikes.txt")
    np.savetxt("lfp.txt", zip(t, lfp), fmt="%18.6f")
