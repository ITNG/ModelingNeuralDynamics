from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_brian_b_jahr_stevens_matches_verified_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    v, B = ns_python.simulate_b_jahr_stevens()
    np.testing.assert_allclose(ns.jahr_stevens_block(v), B, atol=1e-12)
