import time
import numpy as np
import brian2 as b2
from numpy.random import randn
import matplotlib.pyplot as plt
from lib import (plot_data,
                 simulate_PING,
                 tau_peak_function,
                 tau_d_q_function)


if __name__ == "__main__":

    start_time = time.time()

    sigma_e, sigma_i = 0.05, 0.05
    tau_r_e, tau_d_e = 0.5, 3
    tau_r_i, tau_d_i = 0.5, 9
    num_e, num_i = 200, 50
    tau_peak_e = 0.5
    tau_peak_i = 0.5

    g_hat_ee = 0.0  # !
    g_hat_ei = 0.25
    g_hat_ie = 0.25
    g_hat_ii = 0.25

    p_ee = 0.5
    p_ei = 0.5
    p_ie = 0.5
    p_ii = 0.5

    par_e = {
        'num_e': num_e,
        'El': -67.0 * b2.mV,
        'EK': -100.0 * b2.mV,
        'ENa': 50.0 * b2.mV,
        'gl': 0.1 * b2.msiemens,
        'gK': 80.0 * b2.msiemens,
        'gNa': 100.0 * b2.msiemens,
        'C': 1.0 * b2.ufarad,
        'v0': -70 * b2.mV,
        'tau_d_e': tau_d_e * b2.ms,
        'tau_r_e': tau_r_e * b2.ms,
        'tau_peak_e': tau_peak_e * b2.ms
    }

    par_i = {
        'num_i': num_i,
        'El': -65.0 * b2.mV,
        'EK': -90.0 * b2.mV,
        'ENa': 55.0 * b2.mV,
        'gl': 0.1 * b2.msiemens,
        'gK': 9.0 * b2.msiemens,
        'gNa': 35.0 * b2.msiemens,
        'C': 1.0 * b2.ufarad,
        'v0': -63 * b2.mV,
        'tau_d_i': tau_d_i * b2.ms,
        'tau_r_i': tau_r_i * b2.ms,
        'tau_peak_i': tau_peak_i * b2.ms,
    }

    par_syn = {
        'v_rev_e': 0 * b2.mV,
        'v_rev_i': -75 * b2.mV,
        'w_ei': 1,
        'w_ie': 1,
        'w_ee': 1,
        'w_ii': 1,
        'p_ie': p_ie,
        'p_ei': p_ei,
        'p_ee': p_ee,
        'p_ii': p_ii,
    }

    par_e['tau_dq_e'] = tau_d_q_function(tau_d_e,
                                         tau_r_e,
                                         tau_peak_e) * b2.ms
    par_i['tau_dq_i'] = tau_d_q_function(tau_d_i,
                                         tau_r_i,
                                         tau_peak_i) * b2.ms
    par_syn['g_ei'] = g_hat_ei / (num_e * p_ei) * b2.msiemens
    par_syn['g_ee'] = g_hat_ee / (num_e * p_ee) * b2.msiemens
    par_syn['g_ie'] = g_hat_ie / (num_i * p_ei) * b2.msiemens
    par_syn['g_ii'] = g_hat_ii / (num_i * p_ii) * b2.msiemens

    par_e['I_e'] = 1.4 * np.ones(num_e) * (1 + sigma_e * randn(num_e)) * b2.uA
    par_i['I_i'] = 0.0 * np.ones(num_i) * (1 + sigma_i * randn(num_i)) * b2.uA

    par_sim = {
        'simulation_time': 200 * b2.ms,
        'integration_method': 'rk4'
    }
    monitors = simulate_PING(par_e, par_i, par_syn, par_sim)
    print("done in {}".format(time.time()-start_time))
    plot_data(monitors, num_e, num_i)
