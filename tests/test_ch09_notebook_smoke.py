from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter09_notebook_defines_working_examples():
    """Fast (non-MATLAB) sanity check that runs on every push/PR.

    The MATLAB-comparison tests (test_ch09_*.py) are all @pytest.mark.slow,
    so day-to-day CI never parses python/chapter09.ipynb or calls its
    functions -- a broken cell or renamed function would ship green. This
    loads the notebook's definitions and exercises each example with its
    default arguments, checking only shapes and finiteness, not numeric
    accuracy against MATLAB.
    """
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")

    v, w_inf, tau_w = ns.m_current()
    assert v.shape == w_inf.shape == tau_w.shape

    v, r = ns.calcium_rise()
    assert v.shape == r.shape

    t, volt, w = ns.simulate_rtm_m()
    assert t.shape == volt.shape == w.shape
    assert np.all(np.isfinite(volt))

    t, volt, ca = ns.simulate_rtm_ahp()
    assert t.shape == volt.shape == ca.shape
    assert np.all(np.isfinite(volt))

    t, volt, w = ns.simulate_lif_adapt()
    assert t.shape == volt.shape == w.shape
    assert np.all(np.isfinite(volt))

    z, phi = ns.adaptation_map()
    assert z.shape == phi.shape
    assert np.all(np.isfinite(phi))

    t, volt, v_tilde = ns.v_v_tilde()
    assert t.shape == volt.shape == v_tilde.shape
    assert np.all(np.isfinite(volt))
