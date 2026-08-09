from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_lif_voltage_trace_2_matches_analytic_period():
    """Same reasoning as test_ch07_brian_lif_voltage_trace.py: I is tuned
    so the exact period is 20 (time units), verified against the closed
    form rather than MATLAB's non-cumulative plotting arrays.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter07.ipynb")

    tau_m = 2.0
    I = 1 / (1 - np.exp(-20.0 / tau_m)) / tau_m
    t_star = -tau_m * np.log(1 - 1 / (I * tau_m))
    assert abs(t_star - 20.0) < 1e-6

    sm, spm = ns.simulate_LIF_neuron(tau_m=tau_m, I=I, simulation_time=50 * ns.b2.ms)
    spikes = spm.t / ns.b2.ms

    expected = t_star * np.arange(1, len(spikes) + 1)
    assert len(spikes) == 2
    np.testing.assert_allclose(spikes, expected, atol=0.2)
