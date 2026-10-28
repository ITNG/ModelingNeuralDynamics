"""Chapter 39/40 continuous voltage-triggered plasticity in Brian2.

Neither implementation replaces the smooth source dynamics with event rules.
"""

from types import SimpleNamespace
import brian2 as b2
import numpy as np
from .models import make_neurons
from .networks import build_ping


def simulate_short_term(
    *,
    duration=200.0,
    dt=0.01,
    current=0.5,
    facilitation=False,
    U=0.5,
    tau_rec=None,
    tau_facil=500.0,
    U0=0.1,
    u=0.2,
):
    """RTM with dynamic sodium activation and resource/release/gate ODEs.

    p is the available pool, q the active release pool; recovery uses 1-p-q.
    Facilitation evolves log-utilization W, with U=1-exp(-W).
    """
    if not 0 < U < 1 or not 0 < U0 < 1 or not 0 < u < 1:
        raise ValueError("utilizations must lie strictly between zero and one")
    tau_rec = (300.0 if facilitation else 500.0) if tau_rec is None else tau_rec
    equations = """
    dp/dt = (-gamma*p*W+(1-p-q)/tau_rec)/ms : 1
    dq/dt = (gamma*p*W-q/5)/ms : 1
    ds/dt = (q*(1-s)/3-s/9)/ms : 1
    gamma = 1.545*(1+tanh(v/10)) : 1
    """
    equations += (
        "dW/dt = (-(exp(W-W0)-1)/tau_facil+gamma*facil_increment)/ms : 1\n"
        if facilitation
        else "W : 1\n"
    )
    group = make_neurons(
        "rtm",
        current=current,
        dt=dt,
        v0=-70.0,
        instantaneous_m=False,
        extra_equations=equations,
        namespace=dict(
            tau_rec=tau_rec, tau_facil=tau_facil, W0=-np.log1p(-U0), facil_increment=-np.log1p(-u)
        ),
    )
    group.p = 1.0
    group.W = -np.log1p(-(U0 if facilitation else U))
    state = b2.StateMonitor(group, ["v", "p", "q", "s", "W"], record=True)
    spikes = b2.SpikeMonitor(group)
    b2.Network(group, state, spikes).run(duration * b2.ms)
    out = SimpleNamespace(t=np.asarray(state.t / b2.ms), spikes=spikes)
    for name in ("v", "p", "q", "s", "W"):
        setattr(out, name, np.asarray(getattr(state, name)[0]))
    out.U = 1 - np.exp(-out.W)
    return out


def smooth_weight_derivative(v, a, g, *, initial=0.05, C=1.45, tau_plus=10.0, tau_minus=10.0):
    """NumPy diagnostic of the exact chapter 40 row-post/column-pre ODE."""
    g = np.asarray(g)
    mask = g != 0
    kplus = np.where(mask, initial, 0.0)
    kminus = 2 * kplus / 3
    bound = 8 * kplus
    delta = np.maximum(kplus / 2, 1e-6)
    plus = np.exp(-a / tau_plus) * (1 - np.exp(-5 * a / tau_plus))
    minus = np.exp(-a / tau_minus) * (1 - np.exp(-5 * a / tau_minus))
    y = g + kplus * plus[:, None]
    ysoft = np.minimum(bound, y) - delta / 2 * np.log1p(np.exp(-2 * np.abs(bound - y) / delta))
    z = g - kminus * minus[None, :]
    zsoft = np.maximum(0, z) + delta / 2 * np.log1p(np.exp(-2 * np.abs(z) / delta))
    gamma = C * (1 + np.tanh(v / 10))
    return gamma[None, :] * (ysoft - g) + gamma[:, None] * (zsoft - g)


def simulate_three_cell_stdp(*, duration=500.0, dt=0.01, initial=0.05, plastic=True):
    """Two E cells/one I cell with smooth bounded reciprocal E→E STDP.

    Uses the source's currents and chemical conductance strengths. Initial
    voltages have a deterministic tiny heterogeneity instead of phase-splay.
    """
    if initial <= 0:
        raise ValueError("initial weight must be positive")
    r = build_ping(
        n_e=2,
        n_i=1,
        current_e=(0.4, 0.8),
        current_i=0.0,
        g_m=0.0,
        g_ei=0.25,
        g_ie=0.25,
        g_ii=0.25,
        dt=dt,
        adaptation_trace=True,
    )
    r.e.v = -70.0
    r.i.v = -75.0
    # Reset ionic gates to the same initial-voltage steady state.
    for group in (r.e, r.i):
        group.h = "h_inf"
        group.n = "n_inf"
    equations = """
    g_exc_post = weight*s_e_pre : 1 (summed)
    trace_plus = exp(-a_post/10)*(1-exp(-5*a_post/10)) : 1
    trace_minus = exp(-a_pre/10)*(1-exp(-5*a_pre/10)) : 1
    y = weight+kp*trace_plus : 1
    z = weight-km*trace_minus : 1
    soft_upper = (bound+y-abs(bound-y))/2-delta/2*log(1+exp(-2*abs(bound-y)/delta)) : 1
    soft_lower = (z+abs(z))/2+delta/2*log(1+exp(-2*abs(z)/delta)) : 1
    """
    equations += "dweight/dt = learning*1.45*((1+tanh(v_pre/10))*(soft_upper-weight)+(1+tanh(v_post/10))*(soft_lower-weight))/ms : 1\n"
    syn = b2.Synapses(
        r.e,
        r.e,
        equations,
        method="rk4",
        dt=dt * b2.ms,
        namespace=dict(
            kp=initial,
            km=2 * initial / 3,
            bound=8 * initial,
            delta=initial / 2,
            learning=float(plastic),
        ),
    )
    syn.connect(condition="i!=j")
    syn.weight = initial
    r.weight_state = b2.StateMonitor(syn, "weight", record=True)
    r.e_state = b2.StateMonitor(r.e, ["v", "a", "s_e"], record=True)
    r.i_state = b2.StateMonitor(r.i, "v", record=True)
    r.e_spikes = b2.SpikeMonitor(r.e)
    r.i_spikes = b2.SpikeMonitor(r.i)
    r.network.add(syn, r.weight_state, r.e_state, r.i_state, r.e_spikes, r.i_spikes)
    r.network.run(duration * b2.ms)
    r.weights = np.asarray(r.weight_state.weight)
    return r
