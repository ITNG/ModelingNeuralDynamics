from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "22_A_Wilson_Cowan_Model_of_an_Oscillatory_E-I_Network/WILSON_COWAN_LOWERING_W_EE"
MATLAB_DIR = "22/WILSON_COWAN_LOWERING_W_EE"


def test_wilson_cowan_lowering_w_ee_matches_matlab():
    # matlab reuses E/I for all three w_EE runs -- only the last (w_EE=1.0)
    # trace survives to the end of the script
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["E", "I"])

    w_EE, E, I = py.panels[-1]
    assert w_EE == 1.0
    assert np.allclose(E, ref["E"], atol=1e-6)
    assert np.allclose(I, ref["I"], atol=1e-6)
