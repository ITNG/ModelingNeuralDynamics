import brian2 as b2
import matplotlib.pyplot as plt
import numpy as np


def plot_data(monitors, num_e, num_i):

    st_mon_e, st_mon_i, sp_mon_e, sp_mon_i = monitors
    fig, ax = plt.subplots(2, figsize=(10, 5), sharex=True)
    ax[0].plot(sp_mon_e.t / b2.ms, sp_mon_e.i + num_i, '.r', ms=2)
    ax[0].plot(sp_mon_i.t / b2.ms, sp_mon_i.i, '.b', ms=2)

    ax[0].set_xlim(0, np.max(sp_mon_e.t / b2.ms))
    ax[1].set_xlabel('t [ms]', fontsize=14)
    ax[0].set_ylabel('Neuron Index', fontsize=14)
    plt.savefig("PING1.png")
    plt.show()


def simulate_PING(num_e,
                  num_i,
                  I_ext_e,
                  I_ext_i,
                  w_ee,
                  w_ei,
                  w_ie,
                  w_ii,
                  p_ei,
                  p_ie,
                  p_ii,
                  simulation_time,
                  integration_method='euler'):

    # num_rt = 1
    # num_wb = 1
    # I_ext_rt = 1.4 * b2.uA
    # I_ext_wb = 0 * b2.uA

    gSyn_e = 0.25 * b2.msiemens
    gSyn_i = 0.25 * b2.msiemens
    weight_ei = 0.25
    weight_ie = 0.25
    weight_ii = 0.25
    weight_ee = 0.0

    # simulation_time = 200.0 * b2.ms
    # integration_method = "euler"

    # RTM neuron parameters
    El = -67 * b2.mV
    EK = -100 * b2.mV
    ENa = 50 * b2.mV
    ESyn = 0 * b2.mV
    gl = 0.1 * b2.msiemens
    gK = 80 * b2.msiemens
    gNa = 100 * b2.msiemens
    C = 1 * b2.ufarad

    v0_e = -70 * b2.mV
    v0_i = -63 * b2.mV

    tau_r = 0.2 * b2.ms
    tau_d = 2.0 * b2.ms

    # forming RTM model with differential equations
    eqs_e = """

    alphah = 0.128 * exp(-(vm + 50.0*mV) / (18.0*mV))/ms :Hz
    alpham = 0.32/mV * (vm + 54*mV) / (1.0 - exp(-(vm + 54.0*mV) / (4.0*mV)))/ms:Hz
    alphan = 0.032/mV * (vm + 52*mV) / (1.0 - exp(-(vm + 52.0*mV) / (5.0*mV)))/ms:Hz

    betah  = 4.0 / (1.0 + exp(-(vm + 27.0*mV) / (5.0*mV)))/ms:Hz
    betam  = 0.28/mV * (vm + 27.0*mV) / (exp((vm + 27.0*mV) / (5.0*mV)) - 1.0)/ms:Hz
    betan  = 0.5 * exp(-(vm + 57.0*mV) / (40.0*mV))/ms:Hz

    membrane_Im = I_ext + gNa*m**3*h*(ENa-vm) + \
        gl*(El-vm) + gK*n**4*(EK-vm) + gSyn*s_in*(-vm): amp
    I_ext : amp
    s_in  : 1
    gSyn  : siemens

    dm/dt = alpham*(1-m)-betam*m : 1
    dn/dt = alphan*(1-n)-betan*n : 1
    dh/dt = alphah*(1-h)-betah*h : 1
    
    ds/dt = 0.5 * (1 + tanh(0.1*vm/mV)) * (1-s)/tau_r - s/tau_d : 1

    dvm/dt = membrane_Im/C : volt
    """

    neuron_e = b2.NeuronGroup(num_e,
                              eqs_e,
                              method="euler",
                              dt=0.01*b2.ms,
                              threshold='vm>-55*mV',
                              refractory=2*b2.ms)

    # initialize variables
    neuron_e.vm = np.random.rand(num_e) * 10*b2.mV + v0_e
    neuron_e.m = "alpham / (alpham + betam)"
    neuron_e.h = "alphah / (alphah + betah)"
    neuron_e.n = "alphan / (alphan + betan)"
    neuron_e.I_ext = I_ext_e
    neuron_e.s_in = 0
    neuron_e.gSyn = gSyn_e

    # WB neuron parameters
    El = -65 * b2.mV
    EK = -90 * b2.mV
    ENa = 55 * b2.mV
    gl = 0.1 * b2.msiemens
    gK = 9.0 * b2.msiemens
    gNa = 35 * b2.msiemens
    C = 1 * b2.ufarad

    tau_r = 0.5 * b2.ms
    tau_d = 9.0 * b2.ms

    # forming WB model with differential equations
    eqs_i = """
    alphah = 0.35 * exp(-(vm + 58.0*mV) / (20.0*mV))/ms :Hz
    alpham = 0.1/mV * (vm + 35.0*mV) / (1.0 - exp(-0.1/mV * (vm + 35.0*mV))) /ms :Hz
    alphan = -0.05/mV * (vm + 34.0*mV) / (exp(-0.1/mV * (vm + 34.0*mV)) - 1.0)/ms :Hz
    
    betah = 5.0 / (exp(-0.1/mV * (vm + 28.0*mV)) + 1.0)/ms :Hz
    betam = 4.0 * exp(-(vm + 60.0*mV) / (18.0*mV))/ms :Hz
    betan = 0.625 * exp(-(vm + 44.0*mV) / (80.0*mV))/ms :Hz
    
    m_inf = alpham / (alpham + betam) : 1
    
    membrane_Im = I_ext + gNa*m_inf**3*h*(ENa-vm) + \
        gl*(El-vm) + gK*n**4*(EK-vm) +gSyn*s_in*(ESyn-vm): amp
    
    I_ext : amp
    s_in  : 1
    gSyn  : siemens
    
    dn/dt = alphan*(1-n)-betan*n : 1
    dh/dt = alphah*(1-h)-betah*h : 1

    ds/dt = 0.5 * (1 + tanh(0.1*vm/mV)) * (1-s)/tau_r - s/tau_d : 1
    
    dvm/dt = membrane_Im/C : volt
    """
    neuron_i = b2.NeuronGroup(num_i,
                              eqs_i,
                              threshold="vm>-55*mV",
                              refractory=2*b2.ms,
                              method="euler",
                              dt=0.01*b2.ms)

    neuron_i.vm = v0_i
    neuron_i.h = "alphah / (alphah + betah)"
    neuron_i.n = "alphan / (alphan + betan)"
    neuron_i.I_ext = I_ext_i
    neuron_i.s_in = 0
    neuron_i.gSyn = gSyn_i

    # adding Synapses ----------------------------------------------#
    S_ei = b2.Synapses(neuron_e,
                       neuron_i, '''
                    w : 1
                    s_in_post = w*s_pre:1 (summed)
                    ''')
    S_ei.connect('i != j', p=p_ei)
    S_ei.w = weight_ei
    # --------------------------------------------------------------#
    S_ie = b2.Synapses(neuron_i, neuron_e, '''
                    w : 1
                    s_in_post = w*s_pre:1 (summed)
                    ''' )
    S_ie.connect('i != j', p=p_ie)
    S_ie.w = weight_ie
    # --------------------------------------------------------------#
    S_ii = b2.Synapses(neuron_i, neuron_i, '''
                    w : 1
                    s_in_post = w*s_pre:1 (summed)
                    ''' )
    S_ii.connect('i != j', p=p_ii)
    S_ii.w = weight_ii
    # --------------------------------------------------------------#

    # tracking variables
    st_mon_e = b2.StateMonitor(neuron_e, "vm", record=True)
    st_mon_i = b2.StateMonitor(neuron_i, "vm", record=True)

    sp_mon_e = b2.SpikeMonitor(neuron_e)
    sp_mon_i = b2.SpikeMonitor(neuron_i)

    # running the simulation
    net = b2.Network(neuron_e)
    net.add(neuron_i)
    net.add(st_mon_e)
    net.add(st_mon_i)
    net.add(sp_mon_e)
    net.add(sp_mon_i)
    net.add(S_ei)
    net.add(S_ie)
    net.add(S_ii)
    net.run(simulation_time)

    return st_mon_e, st_mon_i, sp_mon_e, sp_mon_i


if __name__ == "__main__":

    num_e = 200
    num_i = 50
    monitors = simulate_PING(num_e=num_e,
                             num_i=num_i,
                             I_ext_e=1.4*b2.uA,
                             I_ext_i=0 * b2.uA,
                             w_ee=0.25,
                             w_ei=0.25,
                             w_ie=0.25,
                             w_ii=0.25,
                             p_ei=0.5,
                             p_ie=0.5,
                             p_ii=0.5,
                             integration_method='euler',
                             simulation_time=100*b2.ms)
    plot_data(monitors, num_e, num_i)
