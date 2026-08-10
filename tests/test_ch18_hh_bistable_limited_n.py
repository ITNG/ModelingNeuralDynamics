from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/HH_BISTABLE_LIMITED_N"
MATLAB_DIR = "18/HH_BISTABLE_LIMITED_N"


@pytest.mark.slow
def test_hh_bistable_limited_n_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(py.v, ref["v"], atol=1e-4)
