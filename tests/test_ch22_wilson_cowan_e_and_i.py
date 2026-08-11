from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_wilson_cowan_e_and_i_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter22.ipynb")
    t, E, I = ns.simulate_wilson_cowan_e_and_i(t_final=50.0)

    assert t.shape == E.shape == I.shape
    assert np.all(np.isfinite(E)) and np.all(np.isfinite(I))
    assert np.all(E >= 0) and np.all(I >= 0)
