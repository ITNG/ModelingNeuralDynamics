"""Continuous voltage-driven PING gates (chapters 32, 33, 38).

These are transparent representative networks, with explicit initial states,
not reproductions of the book's phase-splay population initialization.
"""

from functools import lru_cache
from types import SimpleNamespace
import brian2 as b2
import numpy as np
from .models import make_neurons


@lru_cache(None)
def synaptic_decay(tau_d=3.0, tau_r=0.5, peak=0.5):
    """Calibrate q decay to the source chapter's q=exp(-t/tau_q) gate peak.

    Uses the original midpoint stepping and interpolation at ds/dt=0.
    """

    def peak_time(tau_q):
        s, t, inc = 0.0, 0.0, 1 / tau_r
        for _ in range(100000):
            old_t, old_inc = t, inc
            mid = s + 0.005 * inc
            s += 0.01 * (np.exp(-(t + 0.005) / tau_q) * (1 - mid) / tau_r - mid / tau_d)
            t += 0.01
            inc = np.exp(-t / tau_q) * (1 - s) / tau_r - s / tau_d
            if inc <= 0:
                return (old_t * (-inc) + t * old_inc) / (old_inc - inc)
        raise RuntimeError("gate peak did not converge")

    lo, hi = 1.0, tau_r
    while peak_time(lo) > peak:
        lo /= 2
    while peak_time(hi) < peak:
        hi *= 2
    for _ in range(48):
        mid = (lo + hi) / 2
        if peak_time(mid) <= peak:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def gate_equations(suffix, decay, rise=0.5, peak=0.5):
    tq = synaptic_decay(float(decay), float(rise), float(peak))
    return f"""
    dq_{suffix}/dt = (5*(1+tanh(v/10))*(1-q_{suffix})-q_{suffix}/{tq!r})/ms : 1
    ds_{suffix}/dt = (q_{suffix}*(1-s_{suffix})/{rise!r}-s_{suffix}/{decay!r})/ms : 1
    """


def connect_gate(source, target, gate, destination, conductance, *, condition=None):
    """Each pool writes to its own summed target conductance."""
    syn = b2.Synapses(
        source,
        target,
        f"{destination}_post = weight*{gate}_pre : 1 (summed)\nweight : 1",
        name="continuous_projection*",
    )
    if condition is None:
        syn.connect()
    else:
        syn.connect(condition=condition)
    syn.weight = conductance
    return syn


def build_ping(
    *,
    n_e=40,
    n_i=10,
    g_m=0.5,
    tau_ie=9.0,
    tau_ii=9.0,
    current_e=1.5,
    current_i=0.7,
    g_ei=0.5,
    g_ie=0.5,
    g_ii=0.5,
    dt=0.01,
    seed=63806,
    poisson_rate=0.0,
    poisson_weight=0.02,
    adaptation_trace=False,
    drive_expression=None,
    drive_namespace=None,
):
    """Build an explicit Network; I→E and I→I always use distinct q/s pools.

    Conductance totals are normalized by source count. Independent Poisson
    neurons connect one-to-one to E cells, never broadcasting one event.
    """
    b2.seed(seed)
    rng = np.random.default_rng(seed)
    extra_e = (
        gate_equations("e", 3.0)
        + """
    g_exc : 1
    g_inh : 1
    dw/dt = ((1/(1+exp(-(v+35)/10))-w)*(3.3*exp((v+35)/20)+exp(-(v+35)/20))/400)/ms : 1
    dg_noise/dt = -g_noise/(3*ms) : 1
    """
    )
    if adaptation_trace:
        extra_e += "da/dt=(1-5*a*(1+tanh(v/10)))/ms : 1\n"
    extra_i = gate_equations("ie", tau_ie) + gate_equations("ii", tau_ii) + "g_exc : 1\ng_inh : 1\n"
    e = make_neurons(
        "rtm",
        n_e,
        current=np.broadcast_to(current_e, (n_e,)),
        dt=dt,
        v0=-70 + rng.uniform(-2, 2, n_e),
        extra_equations=extra_e,
        synaptic_current="(g_exc+g_noise)*(0-v)+g_inh*(-75-v)-gm*w*(v-e_k)"
        + ("+" + drive_expression if drive_expression else ""),
        namespace={"gm": g_m, **(drive_namespace or {})},
    )
    i = make_neurons(
        "wb",
        n_i,
        current=np.broadcast_to(current_i, (n_i,)),
        dt=dt,
        v0=-75 + rng.uniform(-2, 2, n_i),
        extra_equations=extra_i,
        synaptic_current="g_exc*(0-v)+g_inh*(-75-v)",
    )
    e.w = "1/(1+exp(-(v+35)/10))"
    if adaptation_trace:
        e.a = 100.0
    projections = [
        connect_gate(e, i, "s_e", "g_exc", g_ei / n_e),
        connect_gate(i, e, "s_ie", "g_inh", g_ie / n_i),
        connect_gate(i, i, "s_ii", "g_inh", g_ii / n_i),
    ]
    objects = [e, i, *projections]
    input_spikes = None
    if poisson_rate:
        noise = b2.PoissonGroup(n_e, rates=poisson_rate * b2.Hz, dt=dt * b2.ms)
        drive = b2.Synapses(
            noise,
            e,
            on_pre="g_noise_post += noise_weight",
            namespace={"noise_weight": poisson_weight},
        )
        drive.connect(j="i")
        input_spikes = b2.SpikeMonitor(noise)
        objects += [noise, drive, input_spikes]
    return SimpleNamespace(e=e, i=i, network=b2.Network(objects), input_spikes=input_spikes)


def simulate_ping(*, duration=200.0, **kwargs):
    r = build_ping(**kwargs)
    evars = ["v", "w", "q_e", "s_e", "g_noise"]
    r.e_state = b2.StateMonitor(r.e, evars, record=True)
    r.i_state = b2.StateMonitor(r.i, ["v", "q_ie", "s_ie", "q_ii", "s_ii"], record=True)
    r.e_spikes = b2.SpikeMonitor(r.e)
    r.i_spikes = b2.SpikeMonitor(r.i)
    r.network.add(r.e_state, r.i_state, r.e_spikes, r.i_spikes)
    r.network.run(duration * b2.ms)
    return r


def simulate_two_cell_coherence(*, duration=200.0, dt=0.01, g_ei=0.35):
    """Chapter 38 GAMMA_COHERENCE_1: 25 ms main/40 ms distractor drives."""
    phase = np.arange(1, 3001) / 3000
    n_main = np.mean(np.expm1(20 * np.cos(np.pi * phase) ** 2))
    n_other = np.mean(np.expm1(0.1 * np.cos(np.pi * phase) ** 2))
    return simulate_ping(
        n_e=1,
        n_i=1,
        g_m=0.0,
        current_e=0.0,
        current_i=0.0,
        g_ei=g_ei,
        g_ie=0.5,
        g_ii=0.15,
        duration=duration,
        dt=dt,
        drive_expression=".5*(exp(20*cos(pi*t/(25*ms))**2)-1)/n_main+.5*(exp(.1*cos(pi*t/(40*ms))**2)-1)/n_other",
        drive_namespace=dict(n_main=n_main, n_other=n_other),
    )
