from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "17/RTM_WITH_M_CURRENT_F_I"


@pytest.mark.slow
def test_rtm_with_m_current_f_i_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    f_forward, f_backward, I_c, I_star, i_ext_vec = ns.simulate_rtm_with_m_current_f_i()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=280)

    # a handful of points sit close to the M-current-induced onset, where
    # tiny cross-language floating-point differences compound over a long
    # (near-threshold) integration into a small but real frequency drift
    assert np.allclose(f_backward, ref["f_vec"], atol=2e-2)
    assert np.isclose(I_c, ref["I_c"], atol=1e-6)
    assert np.isclose(I_star, ref["I_star"], atol=1e-6)
