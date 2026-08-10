from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "31_ING_Rhythms/1_CELL_ING_CONDITION_NUMBERS"


def test_1_cell_ing_condition_numbers_matches_matlab():
    # verified directly against a `matlab -batch make_table` run:
    # percentage_change = -0.7673, -0.3474, -0.4658 (in that order).
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    assert np.isclose(py.pct_i_ext, -0.7673, atol=1e-3)
    assert np.isclose(py.pct_g_ii, -0.3474, atol=1e-3)
    assert np.isclose(py.pct_tau_d, -0.4658, atol=1e-3)
