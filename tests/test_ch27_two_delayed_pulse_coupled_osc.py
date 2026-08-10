from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "27_Phase_Locking_with_Delays/TWO_DELAYED_PULSE_COUPLED_OSC"
MATLAB_DIR = "27/TWO_DELAYED_PULSE_COUPLED_OSC"


@pytest.mark.slow
def test_two_delayed_pulse_coupled_osc_matches_matlab():
    # matlab's script reuses t_spikes_A/t_spikes_B for both delta runs,
    # so by the end of the script they hold only the second (delta=0.7) run.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t_spikes_A", "t_spikes_B"], timeout=30)

    assert np.allclose(py.t_spikes_A_2, ref["t_spikes_A"], atol=1e-9)
    assert np.allclose(py.t_spikes_B_2, ref["t_spikes_B"], atol=1e-9)
