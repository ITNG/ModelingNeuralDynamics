import numpy as np
import brian2 as b2
from numpy import exp
from numpy.random import randn
import matplotlib.pyplot as plt


def plot_data(monitors, num_e, num_i):

    st_mon_e, st_mon_i, sp_mon_e, sp_mon_i = monitors
    fig, ax = plt.subplots(2, figsize=(10, 5), sharex=True)
    ax[0].plot(sp_mon_e.t / b2.ms, sp_mon_e.i + num_i, '.r', ms=2)
    ax[0].plot(sp_mon_i.t / b2.ms, sp_mon_i.i, '.b', ms=2)

    ax[1].plot(st_mon_e.t / b2.ms,
               np.mean(st_mon_e.vm_e / b2.mV, axis=0),
               c='k',
               lw=2)

    ax[0].set_xlim(0, np.max(sp_mon_e.t / b2.ms))
    ax[1].set_xlabel('t [ms]', fontsize=14)
    ax[1].set_ylim(-100, 50)
    ax[0].set_ylabel('Neuron Index', fontsize=14)
    ax[1].set_ylabel('mean(v), E-cells', fontsize=14)
    plt.savefig("PING1.png")
    # plt.show()


def simulate_PING(par_e, par_i, par_syn, par_sim):

    np.random.seed(1)

    num_e = par_e['num_e']
    num_i = par_i['num_i']
    v0_e = par_e['v0']
    v0_i = par_i['v0']

    # forming RTM model with differential equations
    eqs_e = """

    alphah = 0.128 * exp(-(vm_e + 50.0*mV) / (18.0*mV))/ms :Hz
    alpham = 0.32/mV * (vm_e + 54*mV) / (1.0 - exp(-(vm_e + 54.0*mV) / (4.0*mV)))/ms:Hz
    alphan = 0.032/mV * (vm_e + 52*mV) / (1.0 - exp(-(vm_e + 52.0*mV) / (5.0*mV)))/ms:Hz

    betah  = 4.0 / (1.0 + exp(-(vm_e + 27.0*mV) / (5.0*mV)))/ms:Hz
    betam  = 0.28/mV * (vm_e + 27.0*mV) / (exp((vm_e + 27.0*mV) / (5.0*mV)) - 1.0)/ms:Hz
    betan  = 0.5 * exp(-(vm_e + 57.0*mV) / (40.0*mV))/ms:Hz

    membrane_Im = I_ext_e + gNa*m**3*h*(ENa-vm_e) + \
        gl*(El-vm_e) + gK*n**4*(EK-vm_e) + I_syn_e: amp
    I_ext_e   : amp
    I_syn_e_1 : amp
    I_syn_e_2 : amp
    I_syn_e = I_syn_e_1 + I_syn_e_2: amp

    dm/dt = alpham*(1-m)-betam*m : 1
    dn/dt = alphan*(1-n)-betan*n : 1
    dh/dt = alphah*(1-h)-betah*h : 1

    dq/dt = 0.5 * (1+tanh(0.1 * vm_e/mV))*(1-q)*10.0/ms-q/tau_dq_e : 1
    ds_e/dt = q * (1 - s_e)/tau_r_e - s_e/tau_d_e : 1
    
    dvm_e/dt = membrane_Im/C : volt
    """

    # forming WB model with differential equations
    eqs_i = """
    alphah = 0.35 * exp(-(vm_i + 58.0*mV) / (20.0*mV))/ms :Hz
    alpham = 0.1/mV * (vm_i + 35.0*mV) / (1.0 - exp(-0.1/mV * (vm_i + 35.0*mV))) /ms :Hz
    alphan = -0.05/mV * (vm_i + 34.0*mV) / (exp(-0.1/mV * (vm_i + 34.0*mV)) - 1.0)/ms :Hz
    
    betah = 5.0 / (exp(-0.1/mV * (vm_i + 28.0*mV)) + 1.0)/ms :Hz
    betam = 4.0 * exp(-(vm_i + 60.0*mV) / (18.0*mV))/ms :Hz
    betan = 0.625 * exp(-(vm_i + 44.0*mV) / (80.0*mV))/ms :Hz
    
    m_inf = alpham / (alpham + betam) : 1
    
    membrane_Im = I_ext_i + gNa*m_inf**3*h*(ENa-vm_i) + \
        gl*(El-vm_i) + gK*n**4*(EK-vm_i) + I_syn_i: amp
    
    I_ext_i : amp
    I_syn_i_1 : amp
    I_syn_i_2 : amp
    I_syn_i = I_syn_i_1 + I_syn_i_2 : amp

    
    dn/dt = alphan*(1-n)-betan*n : 1
    dh/dt = alphah*(1-h)-betah*h : 1

    #ds_i/dt = 0.5 * (1 + tanh(0.1*vm_i/mV)) * (1-s_i)/tau_r_i - s_i/tau_d_i : 1

    dq/dt = 0.5 * (1+tanh(0.1 * vm_i/mV))*(1-q)*10.0/ms-q/tau_dq_i : 1
    ds_i/dt = q * (1 - s_i)/tau_r_i - s_i/tau_d_i : 1
    
    dvm_i/dt = membrane_Im/C : volt
    """

    neurons_e = b2.NeuronGroup(num_e,
                               eqs_e,
                               method=par_sim['integration_method'],
                               dt=0.05*b2.ms,
                               threshold='vm_e>-20*mV',
                               refractory='vm_e>-20*mV',
                               namespace=par_e)

    # initialize variables
    neurons_e.vm_e = np.random.rand(num_e) * 10*b2.mV + v0_e
    neurons_e.m = "1 / (1 + betam/alpham)"
    neurons_e.h = "1 / (1 + betah/alphah)"
    neurons_e.n = "1 / (1 + betan/alphan)"
    neurons_e.I_ext_e = par_e['I_e']

    #---------------------------------------------------------------#
    neurons_i = b2.NeuronGroup(num_i,
                               eqs_i,
                               method=par_sim['integration_method'],
                               dt=0.05 * b2.ms,
                               threshold='vm_i>-20*mV',
                               refractory='vm_i>-20*mV',
                               namespace=par_i)

    neurons_i.vm_i = np.random.rand(num_i) * 10*b2.mV + v0_i
    neurons_i.h = "1 / (1 + betah / alphah)"
    neurons_i.n = "1 / (1 + betan / alphan)"
    neurons_i.I_ext_i = par_i['I_i']

    # adding Synapses ----------------------------------------------#
    syn_ei_eqs = '''
    w : 1
    I_syn_i_1_post = g_ei * s_e_pre * (v_rev_e - vm_i) :amp (summed)
    '''
    syn_ii_eqs = '''
    w : 1
    I_syn_i_2_post = g_ii * s_i_pre * (v_rev_i - vm_i) :amp (summed)
    '''

    syn_ie_eqs = '''
    w : 1
    I_syn_e_1_post = g_ie * s_i_pre * (v_rev_i - vm_e) :amp (summed)
    '''
    syn_ee_eqs = '''
    w : 1
    I_syn_e_2_post = g_ee * s_e_pre * (v_rev_e - vm_e) :amp (summed)
    '''
    # --------------------------------------------------------------#
    S_ei = b2.Synapses(neurons_e,
                       neurons_i,
                       syn_ei_eqs,
                       namespace=par_syn)
    S_ie = b2.Synapses(neurons_i,
                       neurons_e,
                       syn_ie_eqs,
                       namespace=par_syn)
    S_ii = b2.Synapses(neurons_i,
                       neurons_i,
                       syn_ii_eqs,
                       namespace=par_syn)
    S_ee = b2.Synapses(neurons_e,
                       neurons_e,
                       syn_ee_eqs,
                       namespace=par_syn)

    S_ie.connect('i !=j', p=par_syn['p_ie'])
    S_ie.w = par_syn['w_ie']

    S_ei.connect('i != j', p=par_syn['p_ei'])
    S_ei.w = par_syn['w_ei']

    S_ii.connect('i != j', p=par_syn['p_ii'])
    S_ii.w = par_syn['w_ii']

    S_ee.connect('i != j', p=par_syn['p_ee'])
    S_ee.w = par_syn['w_ee']

    # tracking variables--------------------------------------------#
    st_mon_e = b2.StateMonitor(neurons_e, "vm_e", record=True)
    st_mon_i = b2.StateMonitor(neurons_i, "vm_i", record=True)

    sp_mon_e = b2.SpikeMonitor(neurons_e)
    sp_mon_i = b2.SpikeMonitor(neurons_i)

    # running the simulation ---------------------------------------#
    net = b2.Network(neurons_e)
    net.add(neurons_i)
    net.add(st_mon_e)
    net.add(st_mon_i)
    net.add(sp_mon_e)
    net.add(sp_mon_i)
    net.add(S_ei)
    net.add(S_ie)
    net.add(S_ii)
    net.add(S_ee)
    net.run(par_sim['simulation_time'])

    return st_mon_e, st_mon_i, sp_mon_e, sp_mon_i
# -------------------------------------------------------------------#


def tau_peak_function(tau_d, tau_r, tau_d_q):

    dt = 0.01
    dt05 = 0.5 * dt

    s = 0
    t = 0
    s_inc = exp(-t / tau_d_q) * (1.0 - s) / tau_r - s * tau_d
    while s_inc > 0:
        t_old = t
        s_inc_old = s_inc
        s_tmp = s + dt05 * s_inc
        s_inc_tmp = exp(-(t + dt05) / tau_d_q) * \
            (1.0 - s_tmp) / tau_r - s_tmp / tau_d
        s = s + dt * s_inc_tmp
        t = t + dt
        s_inc = exp(-t / tau_d_q) * (1.0 - s) / tau_r - s / tau_d

    return (t_old * (-s_inc) + t * s_inc_old) / (s_inc_old - s_inc)
# -------------------------------------------------------------------#


def tau_d_q_function(tau_d, tau_r, tau_hat):

    # set an interval for tau_d_q
    tau_d_q_left = 1.0
    while tau_peak_function(tau_d, tau_r, tau_d_q_left) > tau_hat:
        tau_d_q_left *= 0.5

    tau_d_q_right = tau_r
    while tau_peak_function(tau_d, tau_r, tau_d_q_right) < tau_hat:
        tau_d_q_right *= 2.0

    # bisection method
    while tau_d_q_right - tau_d_q_left > 1e-12:
        tau_d_q_mid = 0.5 * (tau_d_q_left + tau_d_q_right)
        if (tau_peak_function(tau_d, tau_r, tau_d_q_mid) <= tau_hat):
            tau_d_q_left = tau_d_q_mid
        else:
            tau_d_q_right = tau_d_q_mid

    return 0.5 * (tau_d_q_left + tau_d_q_right)
# -------------------------------------------------------------------#
