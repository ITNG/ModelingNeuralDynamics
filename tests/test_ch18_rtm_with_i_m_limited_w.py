from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/RTM_WITH_I_M_LIMITED_W"
MATLAB_DIR = "18/RTM_WITH_I_M_LIMITED_W"


def test_rtm_with_i_m_limited_w_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(py.v, ref["v"], atol=1e-3)
