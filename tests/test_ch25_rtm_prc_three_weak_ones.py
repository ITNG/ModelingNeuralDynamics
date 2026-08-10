from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "25_Phase_Response_Curves_(PRCs)/RTM_PRC_THREE_WEAK_ONES"
MATLAB_DIR = "25/RTM_PRC_THREE_WEAK_ONES"


@pytest.mark.slow
def test_rtm_prc_three_weak_ones_matches_matlab():
    # matlab's script ends after the ijk=3 loop iteration, so g_syn/g_vec
    # in the final workspace hold the third (smallest g_syn) run only.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "g_vec", "g_syn", "T"], timeout=60)

    assert np.allclose(py.phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.isclose(py.T, ref["T"], atol=1e-2)
    assert np.isclose(py.g_syn_vec[-1], ref["g_syn"], atol=1e-9)
    assert np.allclose(py.g_vec_list[-1], ref["g_vec"], atol=1e-3)
