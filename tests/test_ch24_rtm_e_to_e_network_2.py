from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "24_Synchronization_by_Fast_Recurrent_Excitation/RTM_E_TO_E_NETWORK_2"
MATLAB_DIR = "24/RTM_E_TO_E_NETWORK_2"


def test_rtm_e_to_e_network_2_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["frequency"], timeout=300)

    assert np.isclose(py.frequency, ref["frequency"], atol=0.5)
