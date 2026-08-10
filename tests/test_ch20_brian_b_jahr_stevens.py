from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]


def test_brian_b_jahr_stevens_matches_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT / "python" / "20_Chemical_Synapses" / "B_JAHR_STEVENS" / "main.py"
    )
    np.testing.assert_allclose(ns.jahr_stevens_block(py.v), py.B, atol=1e-12)
