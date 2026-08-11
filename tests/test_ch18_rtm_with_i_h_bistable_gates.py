from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/RTM_WITH_I_H_BISTABLE_GATES"


@pytest.mark.slow
def test_rtm_with_i_h_bistable_gates_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v, h, n, r, h_star, n_star, r_star, t_final = ns.simulate_rtm_with_i_h_bistable_gates()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["h", "n", "r"], timeout=60)

    assert np.allclose(h, ref["h"], atol=1e-4)
    assert np.allclose(n, ref["n"], atol=1e-4)
    assert np.allclose(r, ref["r"], atol=1e-4)
