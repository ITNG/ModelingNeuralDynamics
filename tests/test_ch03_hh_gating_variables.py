from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter03_notebook_hh_gating_variables():
    """Fast (non-MATLAB) sanity check for python/chapter03.ipynb."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter03.ipynb")

    v, m_inf_v, h_inf_v, n_inf_v, tau_m, tau_h, tau_n = ns.simulate_hh_gating_variables()

    for gate in (m_inf_v, h_inf_v, n_inf_v):
        assert gate.shape == v.shape
        assert np.all(gate >= 0) and np.all(gate <= 1)

    for tau in (tau_m, tau_h, tau_n):
        assert tau.shape == v.shape
        assert np.all(np.isfinite(tau))
        assert np.all(tau > 0)
