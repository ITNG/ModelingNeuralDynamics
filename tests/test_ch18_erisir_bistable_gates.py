from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/ERISIR_BISTABLE_GATES"


@pytest.mark.slow
def test_erisir_bistable_gates_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v, m, h, n, m_star, h_star, n_star, t_final = ns.simulate_erisir_bistable_gates()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["h", "n"])

    assert np.allclose(h, ref["h"], atol=1e-4)
    assert np.allclose(n, ref["n"], atol=1e-4)
