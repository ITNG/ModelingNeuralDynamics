from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "19/INAPIK_SHOW_SLOW_I_K"


@pytest.mark.slow
def test_inapik_show_slow_i_k_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter19.ipynb")
    t, v, n_slow, i_eff = ns.simulate_inapik_show_slow_i_k()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "n_slow"])

    assert np.allclose(v, ref["v"], atol=1e-6)
    assert np.allclose(n_slow, ref["n_slow"], atol=1e-6)
