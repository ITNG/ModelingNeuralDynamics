from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter07_notebook_lif_voltage_trace_2():
    """Fast (non-MATLAB) sanity check for python/chapter07.ipynb (tau_m=2, period-20ms i)."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter07.ipynb")

    tau_m = 2.0
    i = 1 / (1 - np.exp(-20.0 / tau_m)) / tau_m
    t, v = ns.simulate_lif_neuron(tau_m=tau_m, i=i, t_final=50.0)
    assert t.shape == v.shape
    assert np.all(np.isfinite(v))
    assert np.all(v >= 0) and np.all(v <= 1)
