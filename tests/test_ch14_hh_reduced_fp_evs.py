from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "14_Model_Neurons_of_Bifurcation_Type_2/HH_REDUCED_FP_EVS"
MATLAB_DIR = "14/HH_REDUCED_FP_EVS"


def test_hh_reduced_fp_evs_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["real_part", "imag_part"]
    )

    assert np.allclose(py.real_part, ref["real_part"], atol=1e-6)
    assert np.allclose(py.imag_part, ref["imag_part"], atol=1e-6)
