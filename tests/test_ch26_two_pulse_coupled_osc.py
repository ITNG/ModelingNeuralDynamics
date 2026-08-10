from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "26_Phase_Locking_of_Two_Oscillators"
MATLAB_BASE = "26"


def test_two_pulse_coupled_osc_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "TWO_PULSE_COUPLED_OSC" / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "TWO_PULSE_COUPLED_OSC", "make_figure.m",
                             ["t_spikes_A", "t_spikes_B"], timeout=30)

    assert np.allclose(py.t_spikes_A, ref["t_spikes_A"], atol=1e-9)
    assert np.allclose(py.t_spikes_B, ref["t_spikes_B"], atol=1e-9)


def test_two_pulse_coupled_osc_2_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "TWO_PULSE_COUPLED_OSC_2" / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_BASE / "TWO_PULSE_COUPLED_OSC_2", "make_figure.m",
                             ["t_spikes_A", "t_spikes_B"], timeout=30)

    assert np.allclose(py.t_spikes_A, ref["t_spikes_A"], atol=1e-9)
    assert np.allclose(py.t_spikes_B, ref["t_spikes_B"], atol=1e-9)
