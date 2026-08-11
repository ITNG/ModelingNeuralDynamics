from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "17/WB_F_I_CURVE_AT_ONSET"


@pytest.mark.slow
def test_wb_f_i_curve_at_onset_matches_matlab():
    # a couple minutes each side: a few of the 11 points sit close enough
    # to I_c that settling/spiking takes tens of seconds of simulated time
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    f_vec, i_ext_vec, I_c, C, i_ext_low, i_ext_high = ns.simulate_wb_f_i_curve_at_onset()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "C"], timeout=240)

    assert np.allclose(f_vec, ref["f_vec"], atol=1e-2)
    assert np.isclose(I_c, ref["I_c"], atol=1e-4)
    assert np.isclose(C, ref["C"], atol=1e-1)
