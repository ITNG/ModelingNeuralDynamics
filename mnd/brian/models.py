"""Reusable Brian2 conductance neurons in the book's numerical units.

Voltage is numeric mV, current numeric microampere, conductance numeric mS,
and capacitance numeric microfarad. Every time derivative is divided by Brian's
``ms``; ``dt`` is supplied as a numeric number of milliseconds. Explicit clocks
and namespaces keep independently constructed experiments independent.
"""

import brian2 as b2
import numpy as np

MODEL_DEFAULTS = {
    "hh": dict(e_l=-59., e_k=-82., e_na=45., g_l=.3, g_k=36., g_na=120., C=1., n_power=4),
    "rtm": dict(e_l=-67., e_k=-100., e_na=50., g_l=.1, g_k=80., g_na=100., C=1., n_power=4),
    "wb": dict(e_l=-65., e_k=-90., e_na=55., g_l=.1, g_k=9., g_na=35., C=1., n_power=4),
    "erisir": dict(e_l=-70., e_k=-90., e_na=60., g_l=.5, g_k=224., g_na=112., C=1., n_power=2),
}

# Rates are numerical inverse milliseconds. exprel(x)=(exp(x)-1)/x evaluates
# the removable singularities without cancellation or 0/0.
_RATES = {
    "hh": """
    alpha_m = 1/exprel(-(v+45)/10) : 1
    beta_m = 4*exp(-(v+70)/18) : 1
    alpha_h = .07*exp(-(v+70)/20) : 1
    beta_h = 1/(1+exp(-(v+40)/10)) : 1
    alpha_n = .1/exprel(-(v+60)/10) : 1
    beta_n = .125*exp(-(v+70)/80) : 1
    """,
    "rtm": """
    alpha_m = 1.28/exprel(-(v+54)/4) : 1
    beta_m = 1.4/exprel((v+27)/5) : 1
    alpha_h = .128*exp(-(v+50)/18) : 1
    beta_h = 4/(1+exp(-(v+27)/5)) : 1
    alpha_n = .16/exprel(-(v+52)/5) : 1
    beta_n = .5*exp(-(v+57)/40) : 1
    """,
    "wb": """
    alpha_m = 1/exprel(-(v+35)/10) : 1
    beta_m = 4*exp(-(v+60)/18) : 1
    alpha_h = .35*exp(-(v+58)/20) : 1
    beta_h = 5/(1+exp(-(v+28)/10)) : 1
    alpha_n = .5/exprel(-(v+34)/10) : 1
    beta_n = .625*exp(-(v+44)/80) : 1
    """,
    "erisir": """
    alpha_m = 540/exprel((75.5-v)/13.5) : 1
    beta_m = 1.2262*exp(-v/42.248) : 1
    alpha_h = .0035*exp(-v/24.186) : 1
    beta_h = .0884/exprel(-(v+51.25)/5.2) : 1
    alpha_n = 11.8/exprel((95-v)/11.8) : 1
    beta_n = .025*exp(-v/22.222) : 1
    """,
}


def make_neurons(model, n=1, current=0.0, *, dt=0.01, v0=None,
                 instantaneous_m=None, parameters=None, extra_equations="",
                 synaptic_current=None, namespace=None):
    """Construct an unconnected HH, RTM, WB, or Erisir population.

    ``current`` and ``v0`` may be scalars or length-n arrays. Gates start at their
    voltage-dependent equilibria; callers can overwrite them for other source
    initial conditions. ``parameters`` overrides constants from MODEL_DEFAULTS.
    ``extra_equations`` appends adaptation/synapse states; ``synaptic_current``
    optionally supplies an expression for I_syn instead of an assignable state.
    The returned group emits upward-threshold events while preserving continuous
    voltage dynamics (no voltage reset). It must be added to an explicit Network.
    """
    if model not in MODEL_DEFAULTS:
        raise ValueError(f"Unknown neuron model {model!r}; choose {tuple(MODEL_DEFAULTS)}")
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be a positive finite number of milliseconds")
    if not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("n must be a positive integer")
    constants = {**MODEL_DEFAULTS[model], **(parameters or {}), **(namespace or {})}
    if constants["C"] <= 0:
        raise ValueError("C must be positive")
    if instantaneous_m is None:
        instantaneous_m = model != "hh"
    activation = "m_inf" if instantaneous_m else "m"
    sodium_equation = "" if instantaneous_m else "dm/dt = (alpha_m*(1-m)-beta_m*m)/ms : 1"
    current_equation = "I_syn : 1" if synaptic_current is None else f"I_syn = {synaptic_current} : 1"
    equations = _RATES[model] + f"""
    m_inf = alpha_m/(alpha_m+beta_m) : 1
    h_inf = alpha_h/(alpha_h+beta_h) : 1
    n_inf = alpha_n/(alpha_n+beta_n) : 1
    dv/dt = (I_ext + I_syn + g_l*(e_l-v) + g_na*{activation}**3*h*(e_na-v)
             + g_k*n**n_power*(e_k-v))/C/ms : 1
    dh/dt = (alpha_h*(1-h)-beta_h*h)/ms : 1
    dn/dt = (alpha_n*(1-n)-beta_n*n)/ms : 1
    I_ext : 1
    {current_equation}
    {sodium_equation}
    """ + extra_equations
    neuron = b2.NeuronGroup(
        n, equations, method="rk4", threshold="v > -20", refractory="v > -20",
        clock=b2.Clock(dt=dt*b2.ms), namespace=constants, name=f"{model}_neurons*",
    )
    neuron.v = (-63. if model == "wb" else -70.) if v0 is None else v0
    neuron.h = neuron.h_inf[:]
    neuron.n = neuron.n_inf[:]
    if not instantaneous_m:
        neuron.m = neuron.m_inf[:]
    neuron.I_ext = current
    return neuron
