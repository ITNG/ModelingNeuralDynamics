from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "26/F_TILDE"
NOTEBOOK = ROOT / "python" / "chapter26.ipynb"


@pytest.mark.slow
def test_f_tilde_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["s", "phi", "f"], timeout=30)
    s, phi, f = ns.simulate_f_tilde()

    assert np.allclose(s, ref["s"], atol=1e-9)
    assert np.allclose(phi, ref["phi"], atol=1e-9)
    assert np.allclose(f, ref["f"], atol=1e-9)
