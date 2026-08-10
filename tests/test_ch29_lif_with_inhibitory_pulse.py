from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "29_Stability_of_the_Synchronous_State/LIF_WITH_INHIBITORY_PULSE"
MATLAB_DIR = "29/LIF_WITH_INHIBITORY_PULSE"


@pytest.mark.slow
def test_lif_with_inhibitory_pulse_matches_matlab():
    # matlab's script reuses v/period across all N=10 neurons of all 3
    # panels, so by the end of the script they hold only the last
    # neuron's (i=N) trace/period from the third (strong/fast) panel.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["v", "period"], timeout=30)

    v_last, k_last, period_last = py.traces_strong[-1]
    assert np.isclose(period_last, ref["period"], atol=1e-3)
    assert np.allclose(v_last[:k_last + 1], ref["v"][:k_last + 1], atol=1e-3)
