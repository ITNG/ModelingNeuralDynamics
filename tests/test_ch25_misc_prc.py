from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/MISC_PRC"


@pytest.mark.slow
def test_wb_panel_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi_vec, g_vec, T = ns.simulate_wb_prc(i_ext=0.30, g_syn=0.1, v_rev=0.0, N=200, dt=0.01)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR / "WB_PRC", "make_figure.m",
                             ["phi_vec", "g_vec", "T"], timeout=60)

    assert np.allclose(phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.allclose(g_vec, ref["g_vec"], atol=1e-2)


@pytest.mark.slow
def test_hh_panel_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi_vec, g_vec, T = ns.simulate_hh_prc()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR / "HH_PRC", "make_figure.m",
                             ["phi_vec", "g_vec", "T"], timeout=60)

    assert np.allclose(phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.allclose(g_vec, ref["g_vec"], atol=1e-2)


@pytest.mark.slow
def test_erisir_panel_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    phi_vec, g_vec, T = ns.simulate_erisir_prc()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR / "ERISIR_PRC", "make_figure.m",
                             ["phi_vec", "g_vec", "T"], timeout=60)

    assert np.allclose(phi_vec, ref["phi_vec"], atol=1e-9)
    assert np.allclose(g_vec, ref["g_vec"], atol=1e-2)
