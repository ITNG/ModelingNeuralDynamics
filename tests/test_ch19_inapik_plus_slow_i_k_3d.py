from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "19_Bursting/INAPIK_PLUS_SLOW_I_K_3D"
MATLAB_DIR = "19/INAPIK_PLUS_SLOW_I_K_3D"


def test_inapik_plus_slow_i_k_3d_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "n", "n_slow"])

    assert np.allclose(py.v_cyc, ref["v"], atol=1e-6)
    assert np.allclose(py.n_cyc, ref["n"], atol=1e-6)
    assert np.allclose(py.n_slow_cyc, ref["n_slow"], atol=1e-6)
