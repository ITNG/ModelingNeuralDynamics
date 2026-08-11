from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "21/RESET_THRESHOLD"


@pytest.mark.slow
def test_reset_threshold_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter21.ipynb")
    t, v = ns.simulate_reset_threshold()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(v, ref["v"], atol=1e-6)
