from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/ERISIR_BISTABLE"
MATLAB_DIR = "18/ERISIR_BISTABLE"


@pytest.mark.slow
def test_erisir_bistable_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "v_star"])

    assert np.isclose(py.v_star, ref["v_star"], atol=1e-6)
    assert np.allclose(py.v_fire, ref["v"], atol=1e-4)
