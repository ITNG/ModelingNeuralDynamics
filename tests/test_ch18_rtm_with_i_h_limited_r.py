from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/RTM_WITH_I_H_LIMITED_R"


@pytest.mark.slow
def test_rtm_with_i_h_limited_r_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v, t_final, dt = ns.simulate_rtm_with_i_h_limited_r()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"], timeout=60)

    assert np.allclose(v, ref["v"], atol=1e-4)
