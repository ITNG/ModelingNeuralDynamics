from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "24_Synchronization_by_Fast_Recurrent_Excitation/RTM_SYNC"
MATLAB_DIR = "24/RTM_SYNC"


def test_rtm_sync_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t_spikes", "i_spikes"])

    assert np.allclose(py.t_spikes, ref["t_spikes"], atol=1e-6)
    assert np.array_equal(py.i_spikes + 1, ref["i_spikes"].astype(int))
