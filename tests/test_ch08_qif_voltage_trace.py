from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter08_notebook_qif_voltage_trace():
    """Fast (non-MATLAB) sanity check for python/chapter08.ipynb."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter08.ipynb")

    t, v, reset = ns.simulate_qif_voltage_trace()
    assert t.shape == v.shape == reset.shape
    assert np.all(np.isfinite(v))
    assert reset.any()  # the QIF neuron should reset repeatedly at the default drive
