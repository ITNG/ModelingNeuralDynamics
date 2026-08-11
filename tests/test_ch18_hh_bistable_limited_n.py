from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/HH_BISTABLE_LIMITED_N"


@pytest.mark.slow
def test_hh_bistable_limited_n_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v, t_final, dt = ns.simulate_hh_bistable_limited_n()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(v, ref["v"], atol=1e-4)
