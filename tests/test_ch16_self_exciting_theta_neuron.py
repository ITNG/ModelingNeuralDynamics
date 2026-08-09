from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "16_Model_Neurons_of_Bifurcation_Type_3/SELF_EXCITING_THETA_NEURON"
MATLAB_DIR = "16/SELF_EXCITING_THETA_NEURON"


def test_self_exciting_theta_neuron_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["theta", "w"])

    assert np.allclose(py.theta, ref["theta"], atol=1e-6)
    assert np.allclose(py.w, ref["w"], atol=1e-6)
