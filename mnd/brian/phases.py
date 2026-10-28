"""Brian2 protocols for phase resetting, delayed coupling and common inhibition.

Time is numeric milliseconds, except phase-oscillator examples whose intrinsic
period is mapped to one millisecond. Phases are fractions of a natural cycle.
"""
from types import SimpleNamespace

import brian2 as b2
import numpy as np


def _first_spikes(monitor, size):
    first = np.full(size, np.inf)
    np.minimum.at(first, np.asarray(monitor.i), np.asarray(monitor.t/b2.ms))
    if not np.isfinite(first).all():
        raise RuntimeError("Not every oscillator fired within the observation window")
    return first


def theta_prc(phases=None, tau_m=2., current=.13, kick=.1, dt=.002):
    """Measure chapter25's voltage-kick theta PRC, using exact cycle phases.

    The QIF voltage coordinate is v=.5+.5*tan(theta/2); a kick changes this
    coordinate, not theta by a fixed amount. Neurons are initialized from the
    analytical unperturbed orbit, then integrated to their first spike.
    """
    phases = np.linspace(.02, .98, 49) if phases is None else np.asarray(phases, float)
    if np.any((phases <= 0) | (phases >= 1)):
        raise ValueError("phases must lie strictly between zero and one")
    if tau_m*current <= .25:
        raise ValueError("theta PRC requires a regularly firing neuron: tau_m*current > .25")
    scale = np.sqrt(tau_m*current-.25)
    period = np.pi*tau_m/scale
    v0 = .5+scale*np.tan(np.pi*phases-np.pi/2)+kick
    group = b2.NeuronGroup(
        len(phases), "dtheta/dt = (-cos(theta)/tau_m+2*current*(1+cos(theta)))/ms : 1",
        threshold="theta >= pi", reset="theta -= 2*pi", method="rk4",
        clock=b2.Clock(dt=dt*b2.ms), namespace={"tau_m": tau_m, "current": current},
    )
    group.theta = 2*np.arctan(2*v0-1)
    spikes = b2.SpikeMonitor(group)
    b2.Network(group, spikes).run((period+dt)*b2.ms, namespace={})
    first = _first_spikes(spikes, len(phases))
    reset_phase = np.arctan(np.tan(np.pi*phases-np.pi/2)+kick/scale)/np.pi+.5
    return SimpleNamespace(phi=phases, prc=1-phases-first/period,
                           predicted_prc=reset_phase-phases, period=period,
                           first_spikes=first)


def phase_pair(kind="quadratic", epsilon=1., initial_phases=(0., .5),
               delay=0., duration=20., dt=.0005, initial_spike=False):
    """Two source phase oscillators, with actual Brian delayed synaptic events.

    Chapter26 uses g=phi²(1-phi); chapter27 uses epsilon*phi*(1-phi).
    ``initial_spike`` schedules chapter27's A spike at time zero and its delayed
    arrival. Otherwise initial phases are the post-initial-event phases in ch26.
    """
    if kind not in {"quadratic", "delayed", "cubic"}:
        raise ValueError("unknown phase coupling kind")
    if delay < 0 or dt <= 0 or duration <= 0:
        raise ValueError("delay must be nonnegative; dt and duration must be positive")
    formula = {"quadratic": "phi_post**2*(1-phi_post)",
               "delayed": "phi_post*(1-phi_post)",
               "cubic": "2*phi_post*(1-phi_post)**3"}[kind]
    group = b2.NeuronGroup(2, "dphi/dt = 1/ms : 1", threshold="phi >= 1",
                          reset="phi -= 1", method="euler", clock=b2.Clock(dt=dt*b2.ms))
    group.phi = initial_phases
    if initial_spike:
        group.phi[0] = 1.
    coupling = b2.Synapses(group, group, on_pre=f"phi_post += epsilon*({formula})",
                          delay=delay*b2.ms, namespace={"epsilon": epsilon})
    coupling.connect(condition="i != j")
    state = b2.StateMonitor(group, "phi", record=True)
    spikes = b2.SpikeMonitor(group)
    b2.Network(group, coupling, state, spikes).run(duration*b2.ms, namespace={})
    return SimpleNamespace(t=np.asarray(state.t/b2.ms), phase=np.asarray(state.phi),
                           spike_t=np.asarray(spikes.t/b2.ms), spike_i=np.asarray(spikes.i),
                           delay=delay)


def weak_difference(epsilon=.1, initial_phase=.4, detuning=0., duration=None, dt=.01):
    """Integrate chapter28's reduced phase difference ODE with Brian2 RK4."""
    duration = 6/epsilon if duration is None else duration
    if epsilon <= 0:
        raise ValueError("epsilon must be positive")
    group = b2.NeuronGroup(1, """
        p = psi % 1 : 1
        q = (-psi) % 1 : 1
        dpsi/dt = epsilon*(p**2*(1-p)-q**2*(1-q)-detuning)/ms : 1
        """, method="rk4", clock=b2.Clock(dt=dt*b2.ms),
        namespace={"epsilon": epsilon, "detuning": detuning})
    group.psi = initial_phase
    state = b2.StateMonitor(group, "psi", record=True)
    b2.Network(group, state).run(duration*b2.ms, namespace={})
    return SimpleNamespace(t=np.asarray(state.t/b2.ms), psi=np.asarray(state.psi[0]))


def lif_inhibition(strength=.15, tau_syn=9., tau_m=10., current=.12,
                   size=10, duration=60., dt=.001):
    """Chapter29: a common decaying inhibitory conductance and first-spike times.

    No reset is applied: the source experiment observes each cell only until its
    first threshold crossing. Plot callers can truncate traces at first_spikes.
    """
    if current*tau_m <= 1:
        raise ValueError("current*tau_m must exceed one for the baseline oscillator")
    if tau_syn <= 0 or tau_m <= 0:
        raise ValueError("time constants must be positive")
    period = tau_m*np.log(tau_m*current/(tau_m*current-1))
    phases = np.arange(size-1, -1, -1)/size
    group = b2.NeuronGroup(size, """
        dv/dt = (-v/tau_m + current-strength*exp(-(t/ms)/tau_syn)*v)/ms : 1
        """, method="rk4", threshold="v >= 1", refractory="v >= 1",
        clock=b2.Clock(dt=dt*b2.ms),
        namespace={"tau_m": tau_m, "current": current, "strength": strength, "tau_syn": tau_syn})
    group.v = (1-np.exp(-phases*period/tau_m))*tau_m*current
    state = b2.StateMonitor(group, "v", record=True)
    spikes = b2.SpikeMonitor(group)
    b2.Network(group, state, spikes).run(duration*b2.ms, namespace={})
    return SimpleNamespace(t=np.asarray(state.t/b2.ms), v=np.asarray(state.v),
                           first_spikes=_first_spikes(spikes, size), phases=phases, period=period)
