from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "17_Frequency_Current_Curves/SETN_F_I"
MATLAB_DIR = "17/SETN_F_I"


def test_setn_f_i_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["f_low", "f_high"])

    assert np.allclose(py.f_low, ref["f_low"], atol=1e-4)
    assert np.allclose(py.f_high, ref["f_high"], atol=1e-4)
