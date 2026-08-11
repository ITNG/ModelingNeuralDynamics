from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter05_notebook_wb_voltage_trace():
    """Fast (non-MATLAB) sanity check for python/chapter05.ipynb."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter05.ipynb")

    t, v = ns.simulate_wb_voltage_trace()
    assert t.shape == v.shape
    assert np.all(np.isfinite(v))
    assert v.max() > 0  # the WB neuron should spike at the default i_ext=0.75
