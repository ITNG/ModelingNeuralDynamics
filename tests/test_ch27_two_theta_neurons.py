from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "27/TWO_THETA_NEURONS"


@pytest.mark.slow
def test_two_theta_neurons_matches_matlab():
    # very slow both sides (~3.5min in python): 81 grid points x
    # 1e6-step simulations of two delay-coupled theta neurons each.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter27.ipynb")
    epsilon_vec, delta_vec, grid_points = ns.simulate_two_theta_neurons_grid()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["epsilon_vec", "delta_vec"], timeout=1800)

    assert np.allclose(epsilon_vec, ref["epsilon_vec"], atol=1e-9)
    assert np.allclose(delta_vec, ref["delta_vec"], atol=1e-9)

    py_synced = {(round(e, 6), round(d, 6)): s for e, d, s in grid_points}
    # every grid point should agree on sync/desync with the theoretical
    # boundary delta = 1/2 - 1/pi*atan(epsilon/2) (matlab's own criterion
    # for coloring points red vs blue).
    for epsilon, delta, synced in grid_points:
        theoretical = delta > 0.5 - 1 / np.pi * np.arctan(epsilon / 2)
        assert synced == theoretical, (epsilon, delta, synced)
