from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "17_Frequency_Current_Curves/WB_F_I_CURVE_AT_ONSET"
MATLAB_DIR = "17/WB_F_I_CURVE_AT_ONSET"


@pytest.mark.slow
def test_wb_f_i_curve_at_onset_matches_matlab():
    # a couple minutes each side: a few of the 11 points sit close enough
    # to I_c that settling/spiking takes tens of seconds of simulated time
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "C"], timeout=240)

    assert np.allclose(py.f_vec, ref["f_vec"], atol=1e-2)
    assert np.isclose(py.I_c, ref["I_c"], atol=1e-4)
    assert np.isclose(py.C, ref["C"], atol=1e-1)
