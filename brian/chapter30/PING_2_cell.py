'''
code by : Abolfazl Ziaeemehr 
Date    : Aug 2020
Email   : a.ziaeemehr@gmail.com
'''

import brian2 as b2
import matplotlib.pyplot as plt
import numpy as np

# to produce figure 30.3 of Borgers
# using instantaneous rise of synaptic gate variable


def plot_data(st_mon_e, st_mon_i, title=None, c='k'):
    """Plots the state_monitor variables "vm" vs. time.

    Args:
        state_monitor (StateMonitor): the data to plot
        title (string, optional): plot title to display
    """

    fig, ax = plt.subplots(2, figsize=(10, 6), sharex=True)

    ax[0].plot(st_mon_e.t / b2.ms, st_mon_e.vm_e[0] /
               b2.mV, lw=2, c="r", alpha=1, label="neuron RTM")
    ax[0].plot(st_mon_i.t / b2.ms, st_mon_i.vm_i[0] /
               b2.mV, lw=2, c="b", alpha=1, label='neuron WB')

    ax[1].plot(st_mon_e.t / b2.ms, st_mon_e.s_e[0],
               lw=2, c="r", label='s_e')
    ax[1].plot(st_mon_i.t / b2.ms, st_mon_i.s_i[0],
               lw=2, c="b", label='s_i')

    ax[1].set_xlabel("t [ms]", fontsize=14)
    ax[0].set_ylabel("v [mV]", fontsize=14)
    ax[1].set_ylabel("s", fontsize=14)

    ax[0].set_xlim(0, np.max(st_mon_e.t / b2.ms))
    ax[0].set_ylim(-100, 50)
    ax[1].set_ylim(0, 1)
    ax[0].legend()
    ax[1].legend()

    if title is not None:
        ax[0].set_title(title)
    plt.savefig("PING_2_cell.png")
    plt.show()

# Reduced Traub-Miles Model


def simulate_2_cell_PING():

    num_e = 1
    num_i = 1
    I_ext_e = 1.4 * b2.uA
    I_ext_i = 0 * b2.uA

    weight_ei = weight_ie = 1  # connection weight of tr to wb cell.

    v0_e = -70 * b2.mV
    v0_i = -63 * b2.mV

    simulation_time = 200.0 * b2.ms
    integration_method = "rk4"

    params_e = {'El': -67.0 * b2.mV,
                'EK': -100.0 * b2.mV,
                'ENa': 50.0 * b2.mV,
                'gl': 0.1 * b2.msiemens,
                'gK': 80.0 * b2.msiemens,
                'gNa': 100.0 * b2.msiemens,
                'C': 1.0 * b2.ufarad,
                'tau_d_e': 3.0 * b2.ms,
                'tau_r_e': 0.5 * b2.ms,
                # 'tau_dq_e': 0.17228803581986085 * b2.ms
                }

    # WB neuron parameters-------------------------------------------
    params_i = {'El': -65.0 * b2.mV,
                'EK': -90.0 * b2.mV,
                'ENa': 55.0 * b2.mV,
                'gl': 0.1 * b2.msiemens,
                'gK': 9.0 * b2.msiemens,
                'gNa': 35.0 * b2.msiemens,
                'C': 1.0 * b2.ufarad,
                'tau_d_i': 9.0 * b2.ms,
                'tau_r_i': 0.5 * b2.ms,
                # 'tau_dq_i': 0.11628780300799235 * b2.ms,
                }

    params_syn = {
        "g_ei": 0.25 * b2.msiemens,
        "g_ie": 0.25 * b2.msiemens,
        'v_rev_e': 0 * b2.mV,
        'v_rev_i': -75 * b2.mV,
    }

    eqs_e = """

    alphah = 0.128 * exp(-(vm_e + 50.0*mV) / (18.0*mV))/ms :Hz
    alpham = 0.32/mV * (vm_e + 54*mV) / (1.0 - exp(-(vm_e + 54.0*mV) / (4.0*mV)))/ms:Hz
    alphan = 0.032/mV * (vm_e + 52*mV) / (1.0 - exp(-(vm_e + 52.0*mV) / (5.0*mV)))/ms:Hz

    betah  = 4.0 / (1.0 + exp(-(vm_e + 27.0*mV) / (5.0*mV)))/ms:Hz
    betam  = 0.28/mV * (vm_e + 27.0*mV) / (exp((vm_e + 27.0*mV) / (5.0*mV)) - 1.0)/ms:Hz
    betan  = 0.5 * exp(-(vm_e + 57.0*mV) / (40.0*mV))/ms:Hz

    membrane_Im = I_ext + gNa*m**3*h*(ENa-vm_e) + \
        gl*(El-vm_e) + gK*n**4*(EK-vm_e) + I_syn_e: amp
    I_ext : amp
    I_syn_e : amp

    dm/dt = alpham*(1-m)-betam*m : 1
    dn/dt = alphan*(1-n)-betan*n : 1
    dh/dt = alphah*(1-h)-betah*h : 1

    ds_e/dt = 0.5 * (1 + tanh(0.1*vm_e/mV)) * (1-s_e)/tau_r_e - s_e/tau_d_e : 1
    
    dvm_e/dt = membrane_Im/C : volt
    """

    neuron_e = b2.NeuronGroup(num_e,
                              eqs_e,
                              method=integration_method,
                              dt=0.05*b2.ms,
                              threshold='vm_e>-55*mV',
                              refractory='vm_e>-55*mV',
                              namespace=params_e)

    # initialize variables
    neuron_e.vm_e = v0_e
    neuron_e.m = "1 / (1 + betam/alpham)"
    neuron_e.h = "1 / (1 + betah/alphah)"
    neuron_e.n = "1 / (1 + betan/alphan)"
    neuron_e.I_ext = I_ext_e

    # forming WB model with differential equations
    eqs_i = """
    alphah = 0.35 * exp(-(vm_i + 58.0*mV) / (20.0*mV))/ms :Hz
    alpham = 0.1/mV * (vm_i + 35.0*mV) / (1.0 - exp(-0.1/mV * (vm_i + 35.0*mV))) /ms :Hz
    alphan = -0.05/mV * (vm_i + 34.0*mV) / (exp(-0.1/mV * (vm_i + 34.0*mV)) - 1.0)/ms :Hz
    
    betah = 5.0 / (exp(-0.1/mV * (vm_i + 28.0*mV)) + 1.0)/ms :Hz
    betam = 4.0 * exp(-(vm_i + 60.0*mV) / (18.0*mV))/ms :Hz
    betan = 0.625 * exp(-(vm_i + 44.0*mV) / (80.0*mV))/ms :Hz
    
    m_inf = alpham / (alpham + betam) : 1
    
    membrane_Im = I_ext + gNa*m_inf**3*h*(ENa-vm_i) + \
        gl*(El-vm_i) + gK*n**4*(EK-vm_i) + I_syn_i: amp
    
    I_ext : amp
    I_syn_i : amp
    
    dn/dt = alphan*(1-n)-betan*n : 1
    dh/dt = alphah*(1-h)-betah*h : 1

    ds_i/dt = 0.5 * (1 + tanh(0.1*vm_i/mV)) * (1-s_i)/tau_r_i - s_i/tau_d_i : 1
    
    dvm_i/dt = membrane_Im/C : volt
    """

    neuron_i = b2.NeuronGroup(num_i,
                              eqs_i,
                              method=integration_method,
                              dt=0.01 * b2.ms,
                              threshold='vm_i>-55*mV',
                              refractory='vm_i>-55*mV',
                              namespace=params_i)

    neuron_i.vm_i = v0_i
    neuron_i.h = "1 / (1 + betah / alphah)"
    neuron_i.n = "1 / (1 + betan / alphan)"
    neuron_i.I_ext = I_ext_i

    #---------------------------------------------------------------#
    syn_ei_eqs = '''
    w : 1
    I_syn_i_post = g_ei * s_e_pre * (v_rev_e - vm_i) :amp (summed)
    '''
    syn_ie_eqs = '''
    w : 1
    I_syn_e_post = g_ie * s_i_pre * (v_rev_i - vm_e) :amp (summed)
    '''
    #---------------------------------------------------------------#
    S_ei = b2.Synapses(neuron_e,
                       neuron_i,
                       syn_ei_eqs,
                       namespace=params_syn)
    S_ei.connect(i=0, j=0)
    S_ei.w = weight_ei
    #---------------------------------------------------------------#
    S_ie = b2.Synapses(neuron_i,
                       neuron_e,
                       syn_ie_eqs,
                       namespace=params_syn)
    S_ie.connect(i=0, j=0)
    S_ie.w = weight_ie

    # tracking variables
    st_mon_e = b2.StateMonitor(neuron_e, ["vm_e", "s_e"], record=True)
    st_mon_i = b2.StateMonitor(neuron_i, ["vm_i", "s_i"], record=True)

    # running the simulation
    net = b2.Network(neuron_e)
    net.add(neuron_i)
    net.add(st_mon_e)
    net.add(st_mon_i)
    net.add(S_ei)
    net.add(S_ie)
    net.run(simulation_time)
    return st_mon_e, st_mon_i


if __name__ == "__main__":

    st_mon_e, st_mon_wb = simulate_2_cell_PING()
    plot_data(st_mon_e, st_mon_wb)
