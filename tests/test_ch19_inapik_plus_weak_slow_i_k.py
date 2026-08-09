from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "19_Bursting/INAPIK_PLUS_WEAK_SLOW_I_K"
MATLAB_DIR = "19/INAPIK_PLUS_WEAK_SLOW_I_K"


def test_inapik_plus_weak_slow_i_k_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(py.v, ref["v"], atol=1e-6)
