from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "09_Spike_Frequency_Adaptation/M_CURRENT"
MATLAB_DIR = "09/M_CURRENT"


@pytest.mark.slow
def test_m_current_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["w_inf", "tau_w"]
    )

    assert np.allclose(py.w_inf, ref["w_inf"], rtol=1e-6, atol=1e-6)
    assert np.allclose(py.tau_w, ref["tau_w"], rtol=1e-6, atol=1e-6)
