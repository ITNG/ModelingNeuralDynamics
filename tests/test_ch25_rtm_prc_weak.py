from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/RTM_PRC_WEAK"


@pytest.mark.slow
def test_rtm_prc_weak_matches_matlab():
    # matlab's script reuses the variable names phi_vec/g_vec/g_syn for
    # both simulations, so by the end of the script they hold the
    # *second* (g_syn=0.0001) run's results, not the first.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi_vec, g_vec, g_syn_1, phi_vec_2, g_vec_2, g_syn_2, T = ns.simulate_weak_prc_comparison()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "g_vec", "g_syn", "T"], timeout=120)

    assert np.allclose(phi_vec_2, ref["phi_vec"], atol=1e-9)
    assert np.isclose(T, ref["T"], atol=1e-2)
    assert np.isclose(g_syn_2, ref["g_syn"], atol=1e-9)
    assert np.allclose(g_vec_2 / g_syn_2, ref["g_vec"] / ref["g_syn"], atol=1e-1)
