from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter08_notebook_theta_firing():
    """Fast (non-MATLAB) sanity check for python/chapter08.ipynb."""
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter08.ipynb")

    t, theta = ns.simulate_theta_firing()
    assert t.shape == theta.shape
    assert np.all(np.isfinite(theta))
