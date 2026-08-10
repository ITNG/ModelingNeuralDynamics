from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "26_Phase_Locking_of_Two_Oscillators/F_TILDE"
MATLAB_DIR = "26/F_TILDE"


def test_f_tilde_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["s", "phi", "f"], timeout=30)

    assert np.allclose(py.s, ref["s"], atol=1e-9)
    assert np.allclose(py.phi, ref["phi"], atol=1e-9)
    assert np.allclose(py.f, ref["f"], atol=1e-9)
