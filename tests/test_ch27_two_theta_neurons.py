from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "27_Phase_Locking_with_Delays/TWO_THETA_NEURONS"
MATLAB_DIR = "27/TWO_THETA_NEURONS"


@pytest.mark.slow
def test_two_theta_neurons_matches_matlab():
    # very slow both sides (~3.5min in python): 81 grid points x
    # 1e6-step simulations of two delay-coupled theta neurons each.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["epsilon_vec", "delta_vec"], timeout=1800)

    assert np.allclose(py.epsilon_vec, ref["epsilon_vec"], atol=1e-9)
    assert np.allclose(py.delta_vec, ref["delta_vec"], atol=1e-9)

    py_synced = {(round(e, 6), round(d, 6)): s for e, d, s in py.grid_points}
    # every grid point should agree on sync/desync with the theoretical
    # boundary delta = 1/2 - 1/pi*atan(epsilon/2) (matlab's own criterion
    # for coloring points red vs blue).
    for epsilon, delta, synced in py.grid_points:
        theoretical = delta > 0.5 - 1 / np.pi * np.arctan(epsilon / 2)
        assert synced == theoretical, (epsilon, delta, synced)
