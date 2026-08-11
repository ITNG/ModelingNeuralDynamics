from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "19/SQUARE_WAVES"


@pytest.mark.slow
def test_square_waves_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter19.ipynb")
    t, v = ns.simulate_square_waves()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(v, ref["v"], atol=1e-6)
