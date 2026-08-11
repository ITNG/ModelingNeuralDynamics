from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "17/LIF_F_I_CURVE"


@pytest.mark.slow
def test_lif_f_i_curve_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    I, f = ns.simulate_lif_f_i_curve()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["I", "f"])

    assert np.allclose(I, ref["I"], atol=1e-9)
    assert np.allclose(f, ref["f"], atol=1e-6, equal_nan=True)
