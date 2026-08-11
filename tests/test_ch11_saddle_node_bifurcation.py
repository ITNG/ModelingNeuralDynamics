from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "11/SADDLE_NODE_BIFURCATION"


@pytest.mark.slow
def test_fixed_points_match_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter11.ipynb")

    # Closed-form fixed points for panel 1 (a=0.45) -- matlab's x_plus /
    # x_minus (its own subplot(131) doesn't rename these between panels,
    # so this is the a=0.45 solution left in the workspace after the
    # whole script runs... no: the last subplot to compute them is panel 3,
    # but panel 3 has no fixed points and never assigns x_plus/x_minus, so
    # the values left over are panel 2's (a=0.5).
    x_plus, x_minus = ns.saddle_node_fixed_points(a=0.5, b=1.0)

    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["x_plus", "x_minus"]
    )
    assert np.isclose(x_plus, ref["x_plus"], rtol=1e-9)
    assert np.isclose(x_minus, ref["x_minus"], rtol=1e-9)
    # a=b=... at the bifurcation the two fixed points coincide
    assert np.isclose(x_plus, x_minus)
