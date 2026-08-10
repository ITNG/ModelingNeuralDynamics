from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "30_The_PING_Model_of_Gamma_Rhythms/PING_5"


def test_ping_5_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural/statistical properties rather than
    # exact spike times: all three panels should show substantial E and
    # I firing, and panel 1 (dense, all-to-all-ish connectivity) should
    # show a much more regular/periodic population rhythm than panel 2
    # (very sparse connectivity, near-independent firing).
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    assert len(py.t_e_spikes_1) > 500
    assert len(py.t_i_spikes_1) > 100
    assert len(py.t_e_spikes_2) > 500
    assert len(py.t_i_spikes_2) > 100
    assert len(py.t_e_spikes_3) > 500
    assert len(py.t_i_spikes_3) > 100

    def population_rate_cv(t_spikes, bin_width=5.0, t_final=py.t_final):
        bins = np.arange(0, t_final + bin_width, bin_width)
        counts, _ = np.histogram(t_spikes, bins=bins)
        return counts.std() / counts.mean()

    cv_dense = population_rate_cv(py.t_e_spikes_1)
    cv_sparse = population_rate_cv(py.t_e_spikes_2)
    assert cv_dense > cv_sparse
