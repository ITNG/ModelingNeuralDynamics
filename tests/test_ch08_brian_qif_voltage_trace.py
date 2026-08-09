from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_qif_voltage_trace_matches_analytic_period():
    """Same non-cumulative-plotting caveat as the ch07 LIF tests: MATLAB's
    make_figure.m clears its v/t arrays after every spike, so the closed
    form is the only usable reference. QIF_INFINITE_THRESHOLD's own
    Python port already derives it: w=sqrt(tau_m*I-1/4),
    T=2*tau_m/w*atan(1/(2w)) is the time for v to go from 0 to 1, which
    repeats every cycle since each reset returns to v=0.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter08.ipynb")

    tau_m, I = 2.0, 0.15
    w = np.sqrt(tau_m * I - 0.25)
    T = 2 * tau_m / w * np.arctan(1 / (2 * w))

    sm, spm = ns.simulate_QIF_neuron(tau_m=tau_m, I=I, simulation_time=150 * ns.b2.ms)
    spikes = spm.t / ns.b2.ms

    expected = T * np.arange(1, len(spikes) + 1)
    assert len(spikes) == 7
    np.testing.assert_allclose(spikes, expected, atol=0.1)
