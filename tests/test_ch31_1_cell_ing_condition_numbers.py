from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_1_cell_ing_condition_numbers_matches_matlab():
    # verified directly against a `matlab -batch make_table` run:
    # percentage_change = -0.7673, -0.3474, -0.4658 (in that order).
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")
    base_period, pct_i_ext, pct_g_ii, pct_tau_d, v_trace = ns.compute_condition_numbers()

    assert np.isclose(pct_i_ext, -0.7673, atol=1e-3)
    assert np.isclose(pct_g_ii, -0.3474, atol=1e-3)
    assert np.isclose(pct_tau_d, -0.4658, atol=1e-3)
