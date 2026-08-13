from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_BASE = "26"
NOTEBOOK = ROOT / "python" / "chapter26.ipynb"


@pytest.mark.slow
def test_two_pulse_coupled_osc_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "TWO_PULSE_COUPLED_OSC", "make_figure.m",
                             ["t_spikes_A", "t_spikes_B"], timeout=30)
    t_spikes_A, t_spikes_B, _, _ = ns.simulate_two_pulse_coupled_osc(
        ns.g_two_pulse_1, phi_A0=0., phi_B0=0.5)

    assert np.allclose(t_spikes_A, ref["t_spikes_A"], atol=1e-9)
    assert np.allclose(t_spikes_B, ref["t_spikes_B"], atol=1e-9)


@pytest.mark.slow
def test_two_pulse_coupled_osc_2_matches_matlab():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "TWO_PULSE_COUPLED_OSC_2", "make_figure.m",
                             ["t_spikes_A", "t_spikes_B"], timeout=30)
    t_spikes_A, t_spikes_B, _, _ = ns.simulate_two_pulse_coupled_osc(
        ns.g_two_pulse_2, phi_A0=0., phi_B0=0.05)

    assert np.allclose(t_spikes_A, ref["t_spikes_A"], atol=1e-9)
    assert np.allclose(t_spikes_B, ref["t_spikes_B"], atol=1e-9)
