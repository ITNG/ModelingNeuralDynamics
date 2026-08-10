from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "24_Synchronization_by_Fast_Recurrent_Excitation/RTM_TWO_CELL_NETWORK"
MATLAB_DIR = "24/RTM_TWO_CELL_NETWORK"


@pytest.mark.slow
def test_rtm_two_cell_network_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t_spikes", "i_spikes"], timeout=60)

    assert np.array_equal(py.i_spikes, ref["i_spikes"].astype(int))

    # the first 10 spikes (before the deliberate 1e-5 mV perturbation)
    # should match tightly; after that this figure is specifically
    # illustrating an unstable equilibrium, so cross-language floating
    # point differences get exponentially amplified by design -- loosen
    # accordingly (still tight enough to confirm the same qualitative
    # growth-then-plateau behavior, not a divergent trajectory)
    assert np.allclose(py.t_spikes[:10], ref["t_spikes"][:10], atol=1e-3)
    assert np.allclose(py.t_spikes[10:], ref["t_spikes"][10:], atol=1.0)
