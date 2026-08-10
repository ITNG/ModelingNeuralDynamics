from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "20_Chemical_Synapses/RTM_WITH_AUTAPSE_F_I_CURVE"
MATLAB_DIR = "20/RTM_WITH_AUTAPSE_F_I_CURVE"


@pytest.mark.slow
def test_rtm_with_autapse_f_i_curve_matches_matlab():
    # ~4min each side: forward+backward sweep over 31 i_ext values with
    # an autaptic RTM neuron.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=300)

    # near-threshold points (and even the fast-firing tail) accumulate
    # small cross-language floating-point drift over these long near-onset
    # integrations (same effect as ch17/RTM_WITH_M_CURRENT_F_I)
    assert np.allclose(py.f_backward, ref["f_vec"], atol=0.25)
    assert np.isclose(py.I_c, ref["I_c"], atol=1e-6)
    assert np.isclose(py.I_star, ref["I_star"], atol=1e-6)
