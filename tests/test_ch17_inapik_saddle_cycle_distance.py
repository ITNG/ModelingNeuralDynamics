from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "17_Frequency_Current_Curves/INAPIK_SADDLE_CYCLE_DISTANCE"
MATLAB_DIR = "17/INAPIK_SADDLE_CYCLE_DISTANCE"


def test_inapik_saddle_cycle_distance_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["d_vec"])

    assert np.allclose(py.d_vec, ref["d_vec"], atol=1e-6)
