from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "19/INAPIK_PLUS_SLOW_I_K_3D"


@pytest.mark.slow
def test_inapik_plus_slow_i_k_3d_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter19.ipynb")
    v_cyc, n_cyc, n_slow_cyc = ns.simulate_inapik_plus_slow_i_k_3d()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "n", "n_slow"])

    assert np.allclose(v_cyc, ref["v"], atol=1e-6)
    assert np.allclose(n_cyc, ref["n"], atol=1e-6)
    assert np.allclose(n_slow_cyc, ref["n_slow"], atol=1e-6)
