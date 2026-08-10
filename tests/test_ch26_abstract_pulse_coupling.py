from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "26_Phase_Locking_of_Two_Oscillators"
MATLAB_BASE = "26"


@pytest.mark.slow
def test_abstract_pulse_coupling_1_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ABSTRACT_PULSE_COUPLING_1" / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_1", "make_figure.m",
                             ["phi"], timeout=30)

    assert np.allclose(py.phi, ref["phi"], atol=1e-9)
    assert np.allclose(py.g(py.phi), py.f(py.phi) - py.phi, atol=1e-9)
    assert np.allclose(py.bigF(py.phi), py.f(1 - py.phi), atol=1e-9)


@pytest.mark.slow
def test_abstract_pulse_coupling_2_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ABSTRACT_PULSE_COUPLING_2" / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_2", "make_figure.m",
                             ["phi"], timeout=30)

    assert np.allclose(py.phi, ref["phi"], atol=1e-9)
    assert np.allclose(py.g(py.phi), py.epsilon * py.phi * (1 - py.phi) ** 3, atol=1e-9)


@pytest.mark.slow
def test_abstract_pulse_coupling_3_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ABSTRACT_PULSE_COUPLING_3" / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_3", "make_figure.m",
                             ["phi"], timeout=30)

    assert np.allclose(py.phi, ref["phi"], atol=1e-9)


@pytest.mark.slow
def test_abstract_pulse_coupling_4_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ABSTRACT_PULSE_COUPLING_4" / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_4", "make_figure.m",
                             ["phi"], timeout=30)

    assert np.allclose(py.phi, ref["phi"], atol=1e-9)
    assert np.allclose(py.bigG(py.phi), py.phi, atol=1e-6)


@pytest.mark.slow
def test_abstract_pulse_coupling_5_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ABSTRACT_PULSE_COUPLING_5" / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "ABSTRACT_PULSE_COUPLING_5", "make_figure.m",
                             ["phi"], timeout=30)

    assert np.allclose(py.phi, ref["phi"], atol=1e-9)
