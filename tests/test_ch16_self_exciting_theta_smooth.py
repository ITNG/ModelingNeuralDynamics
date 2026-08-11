from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "16/SELF_EXCITING_THETA_SMOOTH"


@pytest.mark.slow
def test_self_exciting_theta_smooth_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter16.ipynb")
    t, theta, w, w_max = ns.simulate_self_exciting_theta_smooth()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["theta", "w"])

    assert np.allclose(theta, ref["theta"], atol=1e-6)
    assert np.allclose(w, ref["w"], atol=1e-6)
