from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "19_Bursting/ERISIR_SHOW_SLOW_I_K"
MATLAB_DIR = "19/ERISIR_SHOW_SLOW_I_K"


def test_erisir_show_slow_i_k_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "n_slow"])

    assert np.allclose(py.v, ref["v"], atol=1e-4)
    assert np.allclose(py.n_slow, ref["n_slow"], atol=1e-4)
