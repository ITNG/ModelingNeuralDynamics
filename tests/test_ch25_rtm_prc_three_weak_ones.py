from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/RTM_PRC_THREE_WEAK_ONES"


@pytest.mark.slow
def test_rtm_prc_three_weak_ones_matches_matlab():
    # matlab's script ends after the ijk=3 loop iteration, so g_syn/g_vec
    # in the final workspace hold the third (smallest g_syn) run only.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi_vec, g_vec_list, g_syn_vec, T = ns.simulate_three_weak_prcs()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi_vec", "g_vec", "g_syn", "T"], timeout=60)

    assert np.allclose(phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.isclose(T, ref["T"], atol=1e-2)
    assert np.isclose(g_syn_vec[-1], ref["g_syn"], atol=1e-9)
    assert np.allclose(g_vec_list[-1], ref["g_vec"], atol=1e-3)
