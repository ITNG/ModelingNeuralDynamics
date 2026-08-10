from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "16_Model_Neurons_of_Bifurcation_Type_3/SETN_PHASE_PLANE"
MATLAB_DIR = "16/SETN_PHASE_PLANE"


@pytest.mark.slow
def test_setn_phase_plane_matches_matlab():
    # matlab overwrites theta/w on every ijk/panel loop -- only the very
    # last trajectory (panel D, ijk=25) survives to the end of the script.
    # Both sides take ~2-3 minutes (dt=1e-4, some trajectories linger for a
    # long simulated time near the near-threshold slow passage).
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["theta", "w"], timeout=240)

    theta, w = py.panels[-1]['trajs'][-1]
    assert np.allclose(theta, ref["theta"], atol=1e-6)
    assert np.allclose(w, ref["w"], atol=1e-6)
