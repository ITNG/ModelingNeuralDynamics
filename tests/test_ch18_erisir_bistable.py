from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/ERISIR_BISTABLE"


@pytest.mark.slow
def test_erisir_bistable_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v_rest, v_fire, v_star, t_final, dt = ns.simulate_erisir_bistable()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "v_star"])

    assert np.isclose(v_star, ref["v_star"], atol=1e-6)
    assert np.allclose(v_fire, ref["v"], atol=1e-4)
