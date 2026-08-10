from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "28_Weakly_Coupled_Oscillators/PLOT_D_TWO_FIXED_POINTS"
MATLAB_DIR = "28/PLOT_D_TWO_FIXED_POINTS"


def test_plot_d_two_fixed_points_matches_matlab():
    # matlab's script reuses the variable name "psi" for the plotting
    # grid ([0:300]/300) and then overwrites it with a scalar in the
    # second bisection loop, so the final workspace "psi" is the stable
    # fixed point, not the grid.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["psi", "psi_0"], timeout=30)

    assert np.isclose(py.psi_0, ref["psi_0"], atol=1e-2)
    assert np.isclose(py.psi_stable, ref["psi"], atol=1e-6)
