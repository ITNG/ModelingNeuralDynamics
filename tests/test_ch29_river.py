from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "29_Stability_of_the_Synchronous_State/RIVER"
MATLAB_DIR = "29/RIVER"


@pytest.mark.slow
def test_river_matches_matlab():
    # matlab's script reuses theta/g_syn for many trajectories, so by the
    # end of the script they hold only the last-computed trajectory
    # (theta(1)=1.25, g_syn(1)=0.25, integrated backwards).
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["theta", "g_syn", "g_ast"], timeout=30)

    assert np.isclose(py.g_ast, ref["g_ast"], atol=1e-3)
    assert np.allclose(py.theta_unstable, ref["theta"], atol=1e-6)
    assert np.allclose(py.g_syn_unstable, ref["g_syn"], atol=1e-6)
