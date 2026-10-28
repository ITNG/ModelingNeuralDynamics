"""Periodic forcing from chapters 35–37, integrated by Brian2.

As in the source chapters, numerical voltage/current values are in mV/uA
and time constants are in ms; LIF voltage is instead normalized to [0, 1].
"""

from types import SimpleNamespace
import brian2 as b2
import numpy as np
from .models import make_neurons


def periodic_shape(t, alpha=5.0, period=25.0, power=2):
    """Mean-one smooth pulse; power=4 gives chapter 37's waveform."""
    phase = np.arange(2000) / 2000
    mean = np.mean(np.expm1(alpha * np.cos(np.pi * phase) ** power))
    return np.expm1(alpha * np.cos(np.pi * np.asarray(t) / period) ** power) / mean


def simulate_periodic_lif(
    current=0.2, *, periodic=True, alpha=5.0, period=25.0, g_bar=0.1, duration=100.0, dt=0.01
):
    """Chapter 35 normalized LIF with periodic or mean tonic inhibition."""
    mean = np.mean(np.expm1(alpha * np.cos(np.pi * np.arange(2000) / 2000) ** 2))
    group = (
        b2.NeuronGroup(
            1,
            """
    dv/dt = (-v/10 + drive - conductance*v)/ms : 1
    conductance = gbar*(exp(alpha*cos(pi*t/(period*ms))**2)-1)/norm : 1
    """,
            threshold="v>1",
            reset="v=0",
            method="rk4",
            dt=dt * b2.ms,
            namespace=dict(
                drive=current,
                gbar=g_bar if periodic else 0.0,
                alpha=alpha,
                period=period,
                norm=mean,
            ),
            name="periodic_lif*",
        )
        if periodic
        else b2.NeuronGroup(
            1,
            "dv/dt=(-v/10+drive-gbar*v)/ms : 1",
            threshold="v>1",
            reset="v=0",
            method="rk4",
            dt=dt * b2.ms,
            namespace=dict(drive=current, gbar=g_bar),
            name="tonic_lif*",
        )
    )
    state, spikes = b2.StateMonitor(group, "v", record=True), b2.SpikeMonitor(group)
    b2.Network(group, state, spikes).run(duration * b2.ms)
    return SimpleNamespace(
        t=np.asarray(state.t / b2.ms),
        v=np.asarray(state.v[0]),
        spike_times=np.asarray(spikes.t / b2.ms),
        rate=len(spikes.t) * 1000 / duration,
    )


def simulate_rtm_forcing(
    currents,
    *,
    kind="pulsed",
    g_bar=0.15,
    alpha=1.0,
    period=25.0,
    power=2,
    duration=1000.0,
    dt=0.01,
):
    """Vectorized independent RTM cells: pulsed current or inhibitory forcing.

    kind is 'constant', 'pulsed' (equal mean current), or 'inhibition'.
    power=4, alpha=1 reproduces the chapter 37 thresholding waveform.
    """
    if kind not in ("constant", "pulsed", "inhibition"):
        raise ValueError("unknown forcing kind")
    currents = np.atleast_1d(np.asarray(currents, dtype=float))
    norm = np.mean(np.expm1(alpha * np.cos(np.pi * np.arange(2000) / 2000) ** power))
    shape = "(exp(alpha*cos(pi*t/(period*ms))**power)-1)/norm"
    expression = (
        "I_ext*(" + shape + "-1)"
        if kind == "pulsed"
        else ("gbar*" + shape + "*(-75-v)" if kind == "inhibition" else "0")
    )
    neurons = make_neurons(
        "rtm",
        len(currents),
        current=currents,
        v0=-70.0,
        dt=dt,
        synaptic_current=expression,
        namespace=dict(alpha=alpha, period=period, power=power, norm=norm, gbar=g_bar),
    )
    state, spikes = b2.StateMonitor(neurons, "v", record=True), b2.SpikeMonitor(neurons)
    b2.Network(neurons, state, spikes).run(duration * b2.ms)
    return SimpleNamespace(
        t=np.asarray(state.t / b2.ms),
        v=np.asarray(state.v),
        currents=currents,
        rates=np.asarray(spikes.count) * 1000 / duration,
        spikes=spikes,
    )
