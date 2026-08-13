from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "28/PLOT_D_TWO_FIXED_POINTS"


@pytest.mark.slow
def test_plot_d_two_fixed_points_matches_matlab():
    # matlab's script reuses the variable name "psi" for the plotting
    # grid ([0:300]/300) and then overwrites it with a scalar in the
    # second bisection loop, so the final workspace "psi" is the stable
    # fixed point, not the grid.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter28.ipynb")
    _, _, _, psi_0, _, psi_stable = ns.find_two_fixed_points()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["psi", "psi_0"], timeout=30)

    assert np.isclose(psi_0, ref["psi_0"], atol=1e-2)
    assert np.isclose(psi_stable, ref["psi"], atol=1e-6)
