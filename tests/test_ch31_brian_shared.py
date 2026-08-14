from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "brian" / "chapter31.ipynb"


def test_ch31_notebook_exposes_shared_helpers():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    py = load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")
    expected = py.tau_d_q_function(9.0, 0.5, 0.5)
    actual = ns.solve_tau_dq(9.0, 0.5, 0.5)
    assert np.isclose(actual, expected, rtol=1e-10)
    assert np.isclose(ns.tau_peak(9.0, 0.5, actual), 0.5, atol=0.02)


def test_probability_validation_is_explicit():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    assert ns.validate_probability("p_ii", 0.5) == 0.5
    with pytest.raises(ValueError, match="p_ii"):
        ns.validate_probability("p_ii", -0.01)
    with pytest.raises(ValueError, match="p_gap"):
        ns.validate_probability("p_gap", 1.01)
