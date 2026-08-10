from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "25_Phase_Response_Curves_(PRCs)/THETA_PRC"
MATLAB_DIR = "25/THETA_PRC"


@pytest.mark.slow
def test_theta_prc_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi", "g"], timeout=30)

    assert np.allclose(py.phi, ref["phi"], atol=1e-9)
    assert np.allclose(py.g, ref["g"], atol=1e-9)
