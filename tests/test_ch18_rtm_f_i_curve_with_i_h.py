from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/RTM_F_I_CURVE_WITH_I_H"
MATLAB_DIR = "18/RTM_F_I_CURVE_WITH_I_H"


def test_rtm_f_i_curve_with_i_h_matches_matlab():
    # very slow both sides (~4-24min): several i_ext points sit close to
    # the two bistability thresholds, where settling/spiking takes a very
    # long simulated time (matlab's while-loop has no cap either).
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=1800)

    assert np.allclose(py.f_backward, ref["f_vec"], atol=2e-2)
    assert np.isclose(py.I_c, ref["I_c"], atol=1e-4)
    assert np.isclose(py.I_star, ref["I_star"], atol=1e-4)
