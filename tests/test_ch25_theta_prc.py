from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/THETA_PRC"


@pytest.mark.slow
def test_theta_prc_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi, g = ns.theta_prc()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi", "g"], timeout=30)

    assert np.allclose(phi, ref["phi"], atol=1e-9)
    assert np.allclose(g, ref["g"], atol=1e-9)
