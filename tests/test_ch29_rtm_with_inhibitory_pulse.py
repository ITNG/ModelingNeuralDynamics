from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "29_Stability_of_the_Synchronous_State/RTM_WITH_INHIBITORY_PULSE"
MATLAB_DIR = "29/RTM_WITH_INHIBITORY_PULSE"


def test_rtm_with_inhibitory_pulse_matches_matlab():
    # matlab's script reuses v/t for all three panels, so by the end of
    # the script they hold only the third (strong/fast) panel's traces.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t", "v"], timeout=60)

    # spiking traces: fixed-step floating-point drift shifts spike
    # timing by a fraction of a step, which blows up a point-by-point
    # comparison right on a spike upstroke even though the traces are
    # otherwise identical; RMSE per neuron is robust to that.
    py_t = ref["t"]  # same grid (same dt/t_final), so no interpolation needed
    for k in range(py.N):
        rmse = trace_rmse(ref["t"], ref["v"][k], py_t, py.v_trace_strong[k])
        assert rmse < 2.0
