from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "17_Frequency_Current_Curves/HH_REDUCED_F_I_CURVE"
MATLAB_DIR = "17/HH_REDUCED_F_I_CURVE"


def test_hh_reduced_f_i_curve_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=180)

    # matlab's f_vec is overwritten in place across all sweeps -- final
    # state is the backward sweep result
    assert np.allclose(py.f_backward, ref["f_vec"], atol=1e-2)
    assert np.isclose(py.I_c, ref["I_c"], atol=1e-6)
    assert np.isclose(py.I_star, ref["I_star"], atol=1e-6)
