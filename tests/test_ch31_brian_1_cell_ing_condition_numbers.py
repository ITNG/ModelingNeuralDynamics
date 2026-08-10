from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]


def test_brian_ing_condition_numbers_match_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_python_port(
        ROOT / "python" / "31_ING_Rhythms" / "1_CELL_ING_CONDITION_NUMBERS" / "main.py"
    )
    result = ns.compute_ing_condition_numbers()
    assert np.isclose(result["pct_i_ext"], py.pct_i_ext, atol=0.08)
    assert np.isclose(result["pct_g_ii"], py.pct_g_ii, atol=0.08)
    assert np.isclose(result["pct_tau_d"], py.pct_tau_d, atol=0.08)
