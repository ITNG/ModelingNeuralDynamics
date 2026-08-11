from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "17/SETN_F_I"


@pytest.mark.slow
def test_setn_f_i_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    f_low, f_high, i_ext_vec = ns.simulate_setn_f_i()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["f_low", "f_high"])

    assert np.allclose(f_low, ref["f_low"], atol=1e-4)
    assert np.allclose(f_high, ref["f_high"], atol=1e-4)
