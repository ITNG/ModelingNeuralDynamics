from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "25_Phase_Response_Curves_(PRCs)/RTM_PRC_SHORT_AND_WEAK"
MATLAB_DIR = "25/RTM_PRC_SHORT_AND_WEAK"


@pytest.mark.slow
def test_rtm_prc_short_and_weak_matches_matlab():
    # matlab's script reuses the variable names g_vec/delta_v for both
    # simulations, so by the end of the script they hold the *second*
    # (delta_v=0.001) run's results, not the first.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "g_vec", "delta_v", "T"], timeout=60)

    assert np.allclose(py.phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.isclose(py.T, ref["T"], atol=1e-2)
    assert np.isclose(py.delta_v_2, ref["delta_v"], atol=1e-9)
    assert np.allclose(py.g_vec_2 / py.delta_v_2, ref["g_vec"] / ref["delta_v"], atol=1e-2)
