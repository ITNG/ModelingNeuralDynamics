from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/RTM_PRC_SHORT"


@pytest.mark.slow
def test_rtm_prc_short_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi_vec, g_vec, T = ns.simulate_voltage_kick_prc(delta_v=4.0)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "g_vec", "T"], timeout=60)

    assert np.allclose(phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.isclose(T, ref["T"], atol=1e-2)
    assert np.allclose(g_vec, ref["g_vec"], atol=1e-3)
