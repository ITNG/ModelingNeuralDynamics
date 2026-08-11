from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/RTM_F_I_CURVE_WITH_I_H"


@pytest.mark.slow
def test_rtm_f_i_curve_with_i_h_matches_matlab():
    # very slow both sides (~4-24min): several i_ext points sit close to
    # the two bistability thresholds, where settling/spiking takes a very
    # long simulated time (matlab's while-loop has no cap either).
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    f_forward, f_backward, I_c, I_star, i_ext_vec = ns.simulate_rtm_f_i_curve_with_i_h()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=1800)

    assert np.allclose(f_backward, ref["f_vec"], atol=2e-2)
    assert np.isclose(I_c, ref["I_c"], atol=1e-4)
    assert np.isclose(I_star, ref["I_star"], atol=1e-4)
