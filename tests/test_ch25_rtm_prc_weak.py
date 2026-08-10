from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "25_Phase_Response_Curves_(PRCs)/RTM_PRC_WEAK"
MATLAB_DIR = "25/RTM_PRC_WEAK"


@pytest.mark.slow
def test_rtm_prc_weak_matches_matlab():
    # matlab's script reuses the variable names phi_vec/g_vec/g_syn for
    # both simulations, so by the end of the script they hold the
    # *second* (g_syn=0.0001) run's results, not the first.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "g_vec", "g_syn", "T"], timeout=120)

    assert np.allclose(py.phi_vec_2, ref["phi_vec"], atol=1e-9)
    assert np.isclose(py.T, ref["T"], atol=1e-2)
    assert np.isclose(py.g_syn_2, ref["g_syn"], atol=1e-9)
    assert np.allclose(py.g_vec_2 / py.g_syn_2, ref["g_vec"] / ref["g_syn"], atol=1e-1)
