from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/THETA_PRC_SHORT_WEAK"


@pytest.mark.slow
def test_theta_prc_short_weak_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi, g_hat = ns.theta_prc_short_weak()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi", "g_hat"], timeout=30)

    assert np.allclose(phi, ref["phi"], atol=1e-9)
    assert np.allclose(g_hat, ref["g_hat"], atol=1e-6)
