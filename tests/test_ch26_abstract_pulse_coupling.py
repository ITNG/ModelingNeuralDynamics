from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_BASE = "26"
NOTEBOOK = ROOT / "python" / "chapter26.ipynb"


@pytest.mark.slow
def test_abstract_pulse_coupling_1_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_1", "make_figure.m",
                             ["phi"], timeout=30)
    phi = ns.abstract_phase_grid()

    assert np.allclose(phi, ref["phi"], atol=1e-9)
    assert np.allclose(ns.g_pulse_1(phi), ns.pulse_map_f(ns.g_pulse_1, phi) - phi, atol=1e-9)
    assert np.allclose(ns.pulse_map_bigF(ns.g_pulse_1, phi), ns.pulse_map_f(ns.g_pulse_1, 1 - phi), atol=1e-9)


@pytest.mark.slow
def test_abstract_pulse_coupling_2_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_2", "make_figure.m",
                             ["phi"], timeout=30)
    phi = ns.abstract_phase_grid()
    epsilon = 2.0

    assert np.allclose(phi, ref["phi"], atol=1e-9)
    assert np.allclose(ns.g_pulse_2(phi, epsilon=epsilon), epsilon * phi * (1 - phi) ** 3, atol=1e-9)


@pytest.mark.slow
def test_abstract_pulse_coupling_3_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_3", "make_figure.m",
                             ["phi"], timeout=30)
    phi = ns.abstract_phase_grid()

    assert np.allclose(phi, ref["phi"], atol=1e-9)


@pytest.mark.slow
def test_abstract_pulse_coupling_4_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_4", "make_figure.m",
                             ["phi"], timeout=30)
    phi = ns.abstract_phase_grid()

    assert np.allclose(phi, ref["phi"], atol=1e-9)
    assert np.allclose(ns.pulse_map_bigG(ns.g_pulse_4, phi), phi, atol=1e-6)


@pytest.mark.slow
def test_abstract_pulse_coupling_5_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_5", "make_figure.m",
                             ["phi"], timeout=30)
    phi = ns.abstract_phase_grid()

    assert np.allclose(phi, ref["phi"], atol=1e-9)
