from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_m_current_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    v, w_inf, tau_w = ns.m_current()

    ref = run_matlab_script(
        ROOT / "matlab" / "09/M_CURRENT", "make_figure.m", ["w_inf", "tau_w"]
    )

    assert np.allclose(w_inf, ref["w_inf"], rtol=1e-6, atol=1e-6)
    assert np.allclose(tau_w, ref["tau_w"], rtol=1e-6, atol=1e-6)
