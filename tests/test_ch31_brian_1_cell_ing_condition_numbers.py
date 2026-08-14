from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_brian_ing_condition_numbers_match_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")
    _, pct_i_ext, pct_g_ii, pct_tau_d, _ = py.compute_condition_numbers()
    result = ns.compute_ing_condition_numbers()
    assert np.isclose(result["pct_i_ext"], pct_i_ext, atol=0.08)
    assert np.isclose(result["pct_g_ii"], pct_g_ii, atol=0.08)
    assert np.isclose(result["pct_tau_d"], pct_tau_d, atol=0.08)
