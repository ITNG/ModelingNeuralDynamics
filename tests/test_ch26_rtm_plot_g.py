from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "26_Phase_Locking_of_Two_Oscillators/RTM_PLOT_G"
MATLAB_DIR = "26/RTM_PLOT_G"


def test_rtm_plot_g_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "bigG_vec", "T"], timeout=60)

    assert np.allclose(py.phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.isclose(py.T, ref["T"], atol=1e-2)
    assert np.allclose(py.bigG_vec, ref["bigG_vec"], atol=1e-2)
