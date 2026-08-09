from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "25_Phase_Response_Curves_(PRCs)/MISC_PRC"
MATLAB_DIR = "25/MISC_PRC"


def test_wb_panel_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR / "WB_PRC", "make_figure.m",
                             ["phi_vec", "g_vec", "T"], timeout=60)

    assert np.allclose(py.phi_vec_wb, ref["phi_vec"], atol=1e-9)
    assert np.allclose(py.g_vec_wb, ref["g_vec"], atol=1e-2)


def test_hh_panel_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR / "HH_PRC", "make_figure.m",
                             ["phi_vec", "g_vec", "T"], timeout=60)

    assert np.allclose(py.phi_vec_hh, ref["phi_vec"], atol=1e-9)
    assert np.allclose(py.g_vec_hh, ref["g_vec"], atol=1e-2)


def test_erisir_panel_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR / "ERISIR_PRC", "make_figure.m",
                             ["phi_vec", "g_vec", "T"], timeout=60)

    assert np.allclose(py.phi_vec_erisir, ref["phi_vec"], atol=1e-9)
    assert np.allclose(py.g_vec_erisir, ref["g_vec"], atol=1e-2)
