from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "17/WB_F_I_CURVE"


@pytest.mark.slow
def test_wb_f_i_curve_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    f_forward, f_backward, I_c, I_star, i_ext_vec = ns.simulate_wb_f_i_curve()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=120)

    assert np.allclose(f_backward, ref["f_vec"], atol=1e-2)
    assert np.isclose(I_c, ref["I_c"], atol=1e-6)
    assert np.isclose(I_star, ref["I_star"], atol=1e-6)
