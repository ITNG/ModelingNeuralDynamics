from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter10_notebook_hh_h_plus_n():
    """Fast (non-MATLAB) sanity check for python/chapter10.ipynb."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter10.ipynb")

    t, n, h = ns.simulate_hh_h_plus_n()
    assert t.shape == n.shape == h.shape
    assert np.all(np.isfinite(n)) and np.all(np.isfinite(h))
