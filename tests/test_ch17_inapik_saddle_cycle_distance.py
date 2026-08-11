from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "17/INAPIK_SADDLE_CYCLE_DISTANCE"


@pytest.mark.slow
def test_inapik_saddle_cycle_distance_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    d_vec, i_ext_vec = ns.simulate_inapik_saddle_cycle_distance()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["d_vec"])

    assert np.allclose(d_vec, ref["d_vec"], atol=1e-6)
