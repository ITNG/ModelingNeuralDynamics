from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_lif_voltage_trace_matches_analytic_period():
    """make_figure.m plots this trace incrementally, clearing its v/t arrays
    after every spike -- so the MATLAB workspace at the end only holds the
    unfinished last segment, not something we can extract spike times from.
    The model is linear and exactly solvable instead: starting every cycle
    at v=0, v(t) = I*tau_m*(1-exp(-t/tau_m)) crosses the threshold v=1 at a
    fixed period t* = -tau_m*ln(1 - 1/(I*tau_m)), repeating identically
    every cycle since each reset returns to the same v=0 starting point.
    """
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter07.ipynb")

    tau_m, I = 10.0, 0.11
    t_star = -tau_m * np.log(1 - 1 / (I * tau_m))

    sm, spm = ns.simulate_LIF_neuron(tau_m=tau_m, I=I, simulation_time=100 * ns.b2.ms)
    spikes = spm.t / ns.b2.ms

    expected = t_star * np.arange(1, len(spikes) + 1)
    assert len(spikes) == 4
    np.testing.assert_allclose(spikes, expected, atol=0.2)
