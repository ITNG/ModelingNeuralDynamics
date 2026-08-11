from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter10_notebook_reduced_hh():
    """Fast (non-MATLAB) sanity check for python/chapter10.ipynb."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter10.ipynb")

    t, v, n, m_inf = ns.simulate_reduced_hh()
    assert t.shape == v.shape == n.shape
    assert np.all(np.isfinite(v))
    assert v.max() > 0  # the reduced HH neuron should spike at the default i_ext=10.0
