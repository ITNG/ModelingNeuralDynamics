from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/HH_BISTABLE_GATES"
MATLAB_DIR = "18/HH_BISTABLE_GATES"


def test_hh_bistable_gates_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["m", "h", "n"])

    assert np.allclose(py.m, ref["m"], atol=1e-4)
    assert np.allclose(py.h, ref["h"], atol=1e-4)
    assert np.allclose(py.n, ref["n"], atol=1e-4)
