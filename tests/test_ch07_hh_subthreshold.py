from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter07_notebook_hh_subthreshold():
    """Fast (non-MATLAB) sanity check for python/chapter07.ipynb.

    Backs LIF_NEURON_WITH_HH, SUBTHR_FOR_HH, and TAU_M_FOR_HH's figures.
    """
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter07.ipynb")

    t, v, m, n, h = ns.simulate_hh_subthreshold()
    assert t.shape == v.shape == m.shape == n.shape == h.shape
    assert np.all(np.isfinite(v))
    for gate in (m, n, h):
        assert np.all(gate >= 0) and np.all(gate <= 1)
