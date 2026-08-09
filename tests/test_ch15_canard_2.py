from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "15_Canard_Explosions/CANARD_2"
MATLAB_DIR = "15/CANARD_2"


def test_canard_2_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["vv", "nn"])

    assert np.allclose(py.vv, ref["vv"], atol=1e-6)
    assert np.allclose(py.nn, ref["nn"], atol=1e-6)
