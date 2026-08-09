from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "17_Frequency_Current_Curves/LIF_F_I_CURVE"
MATLAB_DIR = "17/LIF_F_I_CURVE"


def test_lif_f_i_curve_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["I", "f"])

    assert np.allclose(py.I, ref["I"], atol=1e-9)
    assert np.allclose(py.f, ref["f"], atol=1e-6, equal_nan=True)
