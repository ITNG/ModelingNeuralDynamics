from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "20/B_JAHR_STEVENS"


@pytest.mark.slow
def test_b_jahr_stevens_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    v, B = ns.simulate_b_jahr_stevens()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "B"])

    assert np.allclose(v, ref["v"], atol=1e-9)
    assert np.allclose(B, ref["B"], atol=1e-9)
