from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "09_Spike_Frequency_Adaptation/ADAPTATION_MAP"
MATLAB_DIR = "09/ADAPTATION_MAP"


@pytest.mark.slow
def test_adaptation_map_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["phi"])

    # Same Heun/RK2 scheme and dt as MATLAB (this is an event/threshold-
    # crossing map, not a smooth ODE, so odeint doesn't apply) -- can use a
    # tight-ish tolerance since it's effectively the same discrete algorithm.
    assert np.allclose(py.phi, ref["phi"], rtol=1e-3, atol=1e-4)
