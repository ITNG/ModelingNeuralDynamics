from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "20/RTM_WITH_AUTAPSE_F_I_CURVE"


@pytest.mark.slow
def test_rtm_with_autapse_f_i_curve_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    f_forward, f_backward, I_c, I_star, i_ext_vec = ns.simulate_rtm_with_autapse_f_i_curve()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=300)

    # near-threshold points (and even the fast-firing tail) accumulate
    # small cross-language floating-point drift over these long near-onset
    # integrations (same effect as ch17/RTM_WITH_M_CURRENT_F_I)
    assert np.allclose(f_backward, ref["f_vec"], atol=0.25)
    assert np.isclose(I_c, ref["I_c"], atol=1e-6)
    assert np.isclose(I_star, ref["I_star"], atol=1e-6)
