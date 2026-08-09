from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "17_Frequency_Current_Curves/RTM_WITH_M_CURRENT_F_I"
MATLAB_DIR = "17/RTM_WITH_M_CURRENT_F_I"


def test_rtm_with_m_current_f_i_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=280)

    # a handful of points sit close to the M-current-induced onset, where
    # tiny cross-language floating-point differences compound over a long
    # (near-threshold) integration into a small but real frequency drift
    assert np.allclose(py.f_backward, ref["f_vec"], atol=2e-2)
    assert np.isclose(py.I_c, ref["I_c"], atol=1e-6)
    assert np.isclose(py.I_star, ref["I_star"], atol=1e-6)
