from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "31/ABSTRACT_PULSE_COUPLING_INH_2"


@pytest.mark.slow
def test_abstract_pulse_coupling_inh_2_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")
    varphi = ns.simulate_abstract_pulse_coupling_inh_2()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["varphi"], timeout=30)

    assert np.allclose(varphi, ref["varphi"], atol=1e-9)

    # epsilon/a match the notebook's g_pulse_2 defaults for this example.
    epsilon, a = 0.4, 4.0
    assert np.allclose(ns.g_pulse_2(varphi), -epsilon * varphi * np.tanh((1 - varphi) * a) / np.tanh(a),
                        atol=1e-9)
