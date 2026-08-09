from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "08_Quadratic_Integrate_and_Fire_(QIF)_and_Theta_Neurons/QIF_INFINITE_THRESHOLD"
MATLAB_DIR = "08/QIF_INFINITE_THRESHOLD"


def test_qif_infinite_threshold_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
        ["T", "t_ast", "v_0_to_1", "v_1_to_inf"],
    )

    # Closed-form formulas, not numerical integration -- both sides evaluate
    # the same expression, so this can be tight rather than loose.
    assert np.isclose(py.T, ref["T"], rtol=1e-6)
    assert np.isclose(py.t_ast, ref["t_ast"], rtol=1e-6)
    assert np.allclose(py.v_0_to_1, ref["v_0_to_1"], rtol=1e-6, atol=1e-6)
    assert np.allclose(py.v_1_to_inf, ref["v_1_to_inf"], rtol=1e-6, atol=1e-6)
