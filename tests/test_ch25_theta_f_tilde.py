from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/THETA_F_TILDE"


@pytest.mark.slow
def test_theta_f_tilde_matches_matlab():
    # THETA_F_TILDE is the same underlying computation as THETA_F, with
    # extra plot annotations only -- reuse the same notebook function.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi, f = ns.theta_f()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["phi", "f"])

    assert np.allclose(phi, ref["phi"], atol=1e-9)
    assert np.allclose(f, ref["f"], atol=1e-9)
