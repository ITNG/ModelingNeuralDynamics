from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "27/TWO_DELAYED_PULSE_COUPLED_OSC"


@pytest.mark.slow
def test_two_delayed_pulse_coupled_osc_matches_matlab():
    # matlab's script reuses t_spikes_A/t_spikes_B for both delta runs,
    # so by the end of the script they hold only the second (delta=0.7) run.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter27.ipynb")
    _, _, t_spikes_A_2, t_spikes_B_2 = ns.simulate_two_delayed_pulse_coupled_osc()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t_spikes_A", "t_spikes_B"], timeout=30)

    assert np.allclose(t_spikes_A_2, ref["t_spikes_A"], atol=1e-9)
    assert np.allclose(t_spikes_B_2, ref["t_spikes_B"], atol=1e-9)
