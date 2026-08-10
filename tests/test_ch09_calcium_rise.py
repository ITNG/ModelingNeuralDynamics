from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "09_Spike_Frequency_Adaptation/CALCIUM_RISE"
MATLAB_DIR = "09/CALCIUM_RISE"


@pytest.mark.slow
def test_calcium_rise_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["r"])

    assert np.allclose(py.r, ref["r"], rtol=1e-6, atol=1e-6)
