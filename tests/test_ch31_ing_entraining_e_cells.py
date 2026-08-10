from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "31_ING_Rhythms"


def test_ing_entraining_e_cells_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial E and I
    # firing, most cells active) rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ING_ENTRAINING_E_CELLS" / "main.py")
    assert len(py.t_e_spikes) > 1000
    assert len(py.t_i_spikes) > 200
    counts_e = np.bincount(py.i_e_spikes.astype(int), minlength=py.num_e)
    counts_i = np.bincount(py.i_i_spikes.astype(int), minlength=py.num_i)
    assert np.count_nonzero(counts_e) > 0.9 * py.num_e
    assert np.count_nonzero(counts_i) > 0.9 * py.num_i


def test_ing_entraining_e_cells_2_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ING_ENTRAINING_E_CELLS_2" / "main.py")
    assert len(py.results) == 3
    for t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes in py.results:
        assert len(t_e_spikes) > 500
        assert len(t_i_spikes) > 200
