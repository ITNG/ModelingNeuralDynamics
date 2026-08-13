from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "26/RTM_PLOT_G"
NOTEBOOK = ROOT / "python" / "chapter26.ipynb"


@pytest.mark.slow
def test_rtm_plot_g_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "bigG_vec", "T"], timeout=60)
    phi_vec, bigG_vec, T = ns.simulate_rtm_plot_g()

    assert np.allclose(phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.isclose(T, ref["T"], atol=1e-2)
    assert np.allclose(bigG_vec, ref["bigG_vec"], atol=1e-2)
