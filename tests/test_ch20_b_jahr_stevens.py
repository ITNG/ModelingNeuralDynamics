from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "20_Chemical_Synapses/B_JAHR_STEVENS"
MATLAB_DIR = "20/B_JAHR_STEVENS"


def test_b_jahr_stevens_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "B"])

    assert np.allclose(py.v, ref["v"], atol=1e-9)
    assert np.allclose(py.B, ref["B"], atol=1e-9)
