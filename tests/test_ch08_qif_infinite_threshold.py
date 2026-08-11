from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "08/QIF_INFINITE_THRESHOLD"


@pytest.mark.slow
def test_qif_infinite_threshold_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter08.ipynb")
    T, t_ast, v_0_to_1, v_1_to_inf, v_minus_inf_to_0 = ns.simulate_qif_infinite_threshold()
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
        ["T", "t_ast", "v_0_to_1", "v_1_to_inf"],
    )

    # Closed-form formulas, not numerical integration -- both sides evaluate
    # the same expression, so this can be tight rather than loose.
    assert np.isclose(T, ref["T"], rtol=1e-6)
    assert np.isclose(t_ast, ref["t_ast"], rtol=1e-6)
    assert np.allclose(v_0_to_1, ref["v_0_to_1"], rtol=1e-6, atol=1e-6)
    assert np.allclose(v_1_to_inf, ref["v_1_to_inf"], rtol=1e-6, atol=1e-6)
