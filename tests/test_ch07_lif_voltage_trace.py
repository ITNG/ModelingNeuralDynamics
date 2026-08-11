from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter07_notebook_lif_voltage_trace():
    """Fast (non-MATLAB) sanity check for python/chapter07.ipynb (tau_m=10, i=0.11)."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter07.ipynb")

    t, v = ns.simulate_lif_neuron(tau_m=10.0, i=0.11)
    assert t.shape == v.shape
    assert np.all(np.isfinite(v))
    assert np.all(v >= 0) and np.all(v <= 1)  # resets to 0 once it crosses 1
