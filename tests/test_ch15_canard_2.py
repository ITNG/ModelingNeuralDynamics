from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "15/CANARD_2"


@pytest.mark.slow
def test_canard_2_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter15.ipynb")
    vv, nn, a, I = ns.simulate_canard_2()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["vv", "nn"])

    assert np.allclose(vv, ref["vv"], atol=1e-6)
    assert np.allclose(nn, ref["nn"], atol=1e-6)
