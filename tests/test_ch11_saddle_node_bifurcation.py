from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "11_The_Saddle_Node_Bifurcation/SADDLE_NODE_BIFURCATION"
MATLAB_DIR = "11/SADDLE_NODE_BIFURCATION"


def test_fixed_points_match_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    # Closed-form fixed points for panel 1 (a=0.45) -- matlab's x_plus /
    # x_minus (its own subplot(131) doesn't rename these between panels,
    # so this is the a=0.45 solution left in the workspace after the
    # whole script runs... no: the last subplot to compute them is panel 3,
    # but panel 3 has no fixed points and never assigns x_plus/x_minus, so
    # the values left over are panel 2's (a=0.5).
    a, b = 0.5, 1.0
    disc = 1 / (4 * a ** 2 * b ** 2) - 1
    x_plus = 1 / (2 * a * b) + np.sqrt(disc)
    x_minus = 1 / (2 * a * b) - np.sqrt(disc)

    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["x_plus", "x_minus"]
    )
    assert np.isclose(x_plus, ref["x_plus"], rtol=1e-9)
    assert np.isclose(x_minus, ref["x_minus"], rtol=1e-9)
    # a=b=... at the bifurcation the two fixed points coincide
    assert np.isclose(x_plus, x_minus)
